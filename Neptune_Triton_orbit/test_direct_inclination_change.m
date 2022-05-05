OE_triton = OE_T(end,:)'; 
rv0_triton = rvOrb.orb2rv(OE_triton, const.muN); 
at = OE_triton(1); 
T_triton = 2*pi*sqrt( at^3 / const.muN ); 

OE_aero1 = rvOrb.rv2orb(X_Nsat_e1(end,:), const.muN)'; 
a1 = OE_aero1(1); 
T_aero1 = 2*pi*sqrt( a1^3 / const.muN ); 
rv0_aero1 = X_Nsat_e1(end,:); 

% find common line node for Triton and aerobrake orbit 1 

a2 = const.RN + at/2; 
T_aero2 = 2*pi*sqrt( a2^3 / const.muN ); 
OE_aero2 = OE_aero1; 
OE_aero2(1) = a2; 
rv0_aero2 = rvOrb.orb2rv(OE_aero2, const.muN); 

% propagate orbit 
[t, X_triton] = ode45(@fn.EOM, [0: dt : T_triton], rv0_triton, options); 
[t, X_aero1] = ode45(@fn.EOM, [0: dt : T_aero1], rv0_aero1, options);

% try inclination change - do it AT APOAPSIS. First find apoapsis of
% previous orbit 
i_max_aero1 = rnorm_i_max(X_aero1); 
[t, X_Nsat_ei] = incl_change(X_aero1(i_max_aero1,:), OE_T, dt, const, options); 


% [t, X_Nsat_e2] = incl_change(X_Nsat_ei(end,:), OE_T, dt, const, options); 


% [t, X_Nsat_e3] = incl_change(X_Nsat_e2(end,:), OE_T, dt, const, options); 

% [t, X_aero2] = ode45(@fn.EOM, [0: dt : T_aero2], rv0_aero2, options); 

% plot 
figure()
fn.plot3_xyz(X_triton); 
hold on; grid on; 
fn.plot3_xyz(X_aero1); 
fn.plot3_xyz(X_Nsat_ei); 
% fn.plot3_xyz(X_Nsat_e2); 
% fn.plot3_xyz(X_Nsat_e3); 

legend('Triton', 'aero 1', 'incl 1')

% h = plot3_xyz(X, style, linew, i)

%% subfunctions 

function i_max = rnorm_i_max(X_aero1) 

for i = 1:length(X_aero1) 
    rnorm(i,:) = norm(X_aero1(i, 1:3)); 
end 

i_max = find(rnorm == max(rnorm)); 

end 

function [t, X_Nsat_ei] = incl_change(X_Nsat_ep, OE_T, dt, const, options)

% obtain OE for input state (orbit)
oe_di = rvOrb.rv2orb(X_Nsat_ep(end,:), const.muN); 

% orbit normal 
h = cross(X_Nsat_ep(end, 1:3), X_Nsat_ep(end, 4:6)); 
h = h / norm(h); 

% desired plane change direction 
h_des = h; 

% desired change in inclination = 180 deg - Triton inclination 
di = oe_di(3) - OE_T(end,3); 
if di > 10 * pi/180 
    di = 10 * pi/180; 
elseif di < -10 * pi/180 
    di = -10 * pi/180; 
end 

% line node 
n = cross([0 0 1], h); 
n = n / norm(n); 

% arg of perigee (angle RAAN to periapsis) 
w = oe_di(4); 

% set nu at RAAN (periapsis to RAAN). obtain state at LINE OF NODE 
nu = 2*pi - w; 
oe_di(6) = nu; 
rv_di = rvOrb.orb2rv(oe_di, const.muN); 

% v initial 
vi = rv_di(4:6); 
vi_norm = norm(vi); 

% flight path angle 
e = oe_di(2); 
% nu = oe_di(6); 
nu = oe_di(6) - pi; 
fpa = atan( e*sin(nu) / ( 1 + e*cos(nu) ) ); 

% desired delta velocity magnitude 
dvi_norm = 2 * vi_norm * cos(fpa) * sin(di/2); 
% dvi_norm = 1; 

% desired delta velocity vector 
dvi = h_des * dvi_norm; 

% new state 
oe_di(6) = nu; 
rv_di = rvOrb.orb2rv(oe_di, const.muN); 
rv_di(4:6) = rv_di(4:6) + dvi'; 
oe_di = rvOrb.rv2orb(rv_di, const.muN); 

% fpa again 
fpa = atan( e*sin(nu) / ( 1 + e*cos(nu) ) ); 

% propagate 
T = 2*pi*sqrt( oe_di(1)^3 / const.muN ); 
[t, X_Nsat_ei] = ode45(@fn.EOM, [0: dt : T*0.9], rv_di, options); 

end 



