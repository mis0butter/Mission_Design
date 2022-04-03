function [ell_1, ell_2] = bett_lambert_NT(X_NepSat, t0, phi_des, const, plot_option)

mu_Nep = const.muN; 

abcorr  = 'NONE';

%  Convert the epoch to ephemeris time (secs) 
et_t0   = cspice_str2et( t0 );

% get states --> Neptune to Triton
target   = 'Triton';
frame    = 'ECLIPJ2000';
observer = 'Neptune';
abcorr   = 'NONE';

% get Triton position 
et = et_t0;    % propagate ephemeris time by 1 day in secs 
X_NepTri  = spice_state(et, target, frame, abcorr, observer); 

% get angle between satellite and Triton velocities 
r_dep = X_NepSat(1:3); 
r_arr = X_NepTri(1:3); 
v_S = X_NepSat(4:6); 
v_T = X_NepTri(4:6); 

% angle and velocity angles 
phi_r = acosd( dot(r_dep, r_arr) / (norm(r_dep)*norm(r_arr)) ); 
phi_v = acosd( dot(v_S, v_T) / (norm(v_S)*norm(v_T)) ); 

i = 0; 
X_NT_hist = []; 
while abs(phi_r - phi_des) > 0.01
    
    % propagate days 
    if abs(phi_r - phi_des) > 1
        i = i + 0.01; 
    elseif abs(phi_r - phi_des) > 0.1 
        i = i + 0.001;         
    else 
        i = i + 0.0001; 
    end 
    et = et_t0 + i*86400; 
    
    % get velocity 
    X_NepTri = spice_state(et, target, frame, abcorr, observer); 
    r_arr = X_NepTri(1:3); 
    v_T = X_NepTri(4:6); 
    
    % save 
    X_NT_hist = [X_NT_hist; X_NepTri]; 

    % get angle 
    phi_r = acosd( dot(r_dep, r_arr) / (norm(r_dep)*norm(r_arr)) ); 
    phi_v = acosd( dot(v_S, v_T) / (norm(v_S)*norm(v_T)) ); 
    
end 

% find delta time 
d_et = et - et_t0; 

% units in km 
r_dep = r_dep' ; 
r_arr = r_arr ; 
vd_S = v_S' ; 
va_T = v_T ; 

% propagate to get full Triton orbit 
X_NT_fhist = []; 
target   = 'Triton';
for i = 0 : 0.01 : 6
    
    et = et_t0 + i*86400; 
    
    % get state 
    X_NepTri  = spice_state(et, target, frame, abcorr, observer); 
    
    % save Mars vector 
    X_NT_fhist = [X_NT_fhist; X_NepTri]; 
end 

OE0 = rvOrb.rv2orb(X_NepSat, const.muN); 
T0 = 2*pi*sqrt( OE0(1)^3 / const.muN ); 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% propagate to get full (initial) satellite orbit 
[t, X_NS_fhist] = ode45(@fn.EOM, [0 T0], X_NepSat, options); 

% propagate to get how far the satellite would have traveled, no burn 
[t, X_NS_hist] = ode45(@fn.EOM, [0 d_et], X_NepSat, options); 

%% bettadpur lambert 

% Lambert - from VALLADO ed. 4, pgs 467 - 475 

% Departure (Neptune-satellite)
r_dep_mag = norm(r_dep); 

% Arrival (Neptune-Triton) 
r_arr_mag = norm(r_arr); 

% Vallado method ... 
cos_dv = dot(r_dep, r_arr) / (r_dep_mag * r_arr_mag); 

% chord 
c = sqrt( r_dep_mag^2 + r_arr_mag^2 - 2*r_dep_mag*r_arr_mag*cos_dv ); 

% semiperimeter 
s = ( r_arr_mag + r_dep_mag + c ) / 2; 

% min semimajor axis 
amin = s/2;  

[ell_1, ell_2] = a2tof(amin, mu_Nep, s, c, r_dep, r_arr); 

% % mean motion 
% n  = sqrt(mu_Nep/amin^3); 
% 
% % time of flight 
% alpha1 = 2 * asin( sqrt( s/(2*amin) ) ); 
% alpha2 = 2*pi - alpha1; 
% beta1  = 2 * asin( sqrt( (s-c)/(2*amin) ) ); 
% beta2  = - beta1; 
% 
% dE1 = alpha1 - beta1; 
% dE2 = alpha2 - beta2; 
% 
% % solution 4 
% t_a1 = alpha1 - sin(alpha1); 
% t_a2 = alpha2 - sin(alpha2); 
% t_b1 = beta1 - sin(beta1); 
% t_b2 = beta2 - sin(beta2); 
% 
% d_t1 = 1/n * (t_a1 + t_b1); 
% d_t2 = 1/n * (t_a2 + t_b2); 
% 
% f1 = 1 - amin / r_dep_mag * ( 1 - cos(dE1) ); 
% f2 = 1 - amin / r_dep_mag * ( 1 - cos(dE2) ); 
% g1 = d_t1 - 1/n * ( dE1 - sin(dE1) );  
% g2 = d_t2 - 1/n * ( dE2 - sin(dE2) );  
% 
% v_E1 = 1/g1 * ( r_NepTri - f1 * X_NepSat );
% v_E2 = 1/g2 * ( r_NepTri - f2 * X_NepSat );
% 
% % ellipse 1 - rv is from Sun to satellite 
% rv0_e1 = [X_NepSat, v_E1]; 
% rv0_e2 = [X_NepSat, v_E2]; 
% 
% ell_1.rv0 = rv0_e1; 
% ell_1.tof = d_t1; 
% ell_1.amin = amin; 
% 
% ell_2.rv0 = rv0_e2; 
% ell_2.tof = d_t2; 
% ell_2.amin = amin; 
% 
% oe_e1 = rvOrb.rv2orb(rv0_e1', mu_Nep); 
% oe_e2 = rvOrb.rv2orb(rv0_e2', mu_Nep); 

%% propagate 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% ellipse 1 
X_start = ell_1.rv0; 
t_end = ell_1.tof; 
% propagate 
[t, X_hist] = ode45(@fn.EOM, [0 t_end], X_start, options); 
X_ell_1 = X_hist; 

% ellipse 2 
X_start = ell_2.rv0; 
t_end = ell_2.tof; 
% propagate 
[t, X_hist] = ode45(@fn.EOM, [0 t_end], X_start, options); 
X_ell_2 = X_hist; 

%% plot 

if plot_option > 0 
    
    figure()
    plot3(X_ell_1(:,1), X_ell_1(:,2), X_ell_1(:,3), 'b'); 
    hold on; grid on; 
    plot3(X_ell_2(:,1), X_ell_2(:,2), X_ell_2(:,3), 'r--'); 
    
end 

end 

function [ell_1, ell_2] = a2tof(a, mu_Nep, s, c, r_dep, r_arr)

% Departure (Neptune-satellite)
r_dep_mag = norm(r_dep); 

% Arrival (Neptune-Triton) 
r_arr_mag = norm(r_arr); 

% mean motion 
n  = sqrt(mu_Nep/a^3); 

% time of flight 
alpha1 = 2 * asin( sqrt( s/(2*a) ) ); 
alpha2 = 2*pi - alpha1; 
beta1  = 2 * asin( sqrt( (s-c)/(2*a) ) ); 
beta2  = - beta1; 

dE1 = alpha1 - beta1; 
dE2 = alpha2 - beta2; 

% solution 4 
t_a1 = alpha1 - sin(alpha1); 
t_a2 = alpha2 - sin(alpha2); 
t_b1 = beta1 - sin(beta1); 
t_b2 = beta2 - sin(beta2); 

% period 
T = 2*pi*sqrt( a^3 / mu_Nep ); 

d_t1 = 1/n * (t_a1 + t_b1); 
d_t2 = T - d_t1; 

f1 = 1 - a / r_dep_mag * ( 1 - cos(dE1) ); 
f2 = 1 - a / r_dep_mag * ( 1 - cos(dE2) ); 
g1 = d_t1 - 1/n * ( dE1 - sin(dE1) );  
g2 = d_t2 - 1/n * ( dE2 - sin(dE2) );  

v_E1 = 1/g1 * ( r_arr - f1 * r_dep );
v_E2 = 1/g2 * ( r_arr - f2 * r_dep );

% ellipse 1 - rv is from Sun to satellite 
rv0_e1 = [r_dep, v_E1]; 
rv0_e2 = [r_dep, v_E2]; 

ell_1.rv0 = rv0_e1; 
ell_1.tof = d_t1; 
ell_1.amin = a; 

ell_2.rv0 = rv0_e2; 
ell_2.tof = d_t2; 
ell_2.amin = a; 

oe_e1 = rvOrb.rv2orb(rv0_e1', mu_Nep); 
oe_e2 = rvOrb.rv2orb(rv0_e2', mu_Nep); 

end 