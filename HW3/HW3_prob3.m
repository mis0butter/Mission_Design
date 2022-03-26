%% HW 3
% Junette Hsin 

% Hohmann transfer - Vallado 

% Sphere of influence - Vallado / Russell's notes 

% close all; 
clear; 
addpath(genpath('Lambert Battin'))
addpath(genpath('../mice')); 
addpath(genpath('spice_data')); 

%  Load kernel file 
cspice_furnsh( 'spice_data/naif0011.tls' )
cspice_furnsh( 'spice_data/de421.bsp' )       
cspice_furnsh( 'spice_data/pck00010.tpc ') 

% sun mu (m^3/s^2)
mu_sun_m = 1.32712440018e20; 
mu_sun_km = mu_sun_m / (1000^3); 
mu = mu_sun_km; 

%% Hohmann Transfer 

%  Define parameters for a state lookup:
% t0      = 'Oct 20, 2020 11:00 AM CST'; 
t0      = 'May 22, 2020'; 

%  Convert the epoch to ephemeris time (secs) 
et_t0   = cspice_str2et( t0 );

abcorr  = 'NONE';

% get states --> Sun to Earth 
target   = 'Earth';
frame    = 'J2000'; 
observer = 'Sun';
abcorr   = 'NONE'; 
et = et_t0;    

% get Sun-earth position 
X_sunE  = spice_state(et, target, frame, abcorr, observer); 

% get sun-Mars position 
target   = 'Mars';
X_sunM  = spice_state(et, target, frame, abcorr, observer); 

% get angle between Earth and Mars velocities 
r_sunE = X_sunE(1:3); 
r_sunM = X_sunM(1:3); 
v_sunE = X_sunE(4:6); 
v_sunM = X_sunM(4:6); 

% position and velocity angles 
phi_r = acosd( dot(r_sunE, r_sunM) / (norm(r_sunE)*norm(r_sunM)) ); 
phi_v = acosd( dot(v_sunE, v_sunM) / (norm(v_sunE)*norm(v_sunM)) ); 

% desired phase angle 
phi_des = 179.6234; 

i = 0; 
while abs(phi_r - phi_des) > 0.01
    
    if abs(phi_r - phi_des) > 10 
        % propagate by 1 day 
        i = i + 1;         
    elseif abs(phi_r - phi_des) > 1 
        % propagate by 0.1 day 
        i = i + 0.1; 
    else
        % propagate by 0.01 day 
        i = i + 0.01;         
    end 
    et = et_t0 + i*86400; 
    
    % get velocity 
    X_sunM  = spice_state(et, target, frame, abcorr, observer); 
    r_sunM = X_sunM(1:3); 
    v_sunM = X_sunM(4:6); 

    % get angle 
    phi_r = acosd( dot(r_sunE, r_sunM) / (norm(r_sunE)*norm(r_sunM)) )
    phi_v = acosd( dot(v_sunE, v_sunM) / (norm(v_sunE)*norm(v_sunM)) ); 
    
end 

%% bettadpur lambert 

% Lambert - from VALLADO ed. 4, pgs 467 - 475 

% Arrival (Mars), AU units 
rM_mag = norm(r_sunM); 
rE_mag = norm(r_sunE); 

% Vallado method ... 
cos_dv = dot(r_sunE, r_sunM) / (rE_mag * rM_mag); 

% chord 
c = sqrt( rE_mag^2 + rM_mag^2 - 2*rE_mag*rM_mag*cos_dv ); 

% semiperimeter 
s = ( rM_mag + rE_mag + c ) / 2; 

% min semimajor axis 
amin = s/2;  

% mean motion 
n  = sqrt(mu/amin^3); 

% time of flight 
alpha1 = 2 * asin( sqrt( s/(2*amin) ) ); 
alpha2 = 2*pi - alpha1; 
beta1  = 2 * asin( sqrt( (s-c)/(2*amin) ) ); 
beta2  = - beta1; 

dE1 = alpha1 - beta1; 
dE2 = alpha2 - beta2; 

% solution 4 
t_a1 = alpha1 - sin(alpha1); 
t_a2 = alpha2 - sin(alpha2); 
t_b1 = beta1 - sin(beta1); 
t_b2 = beta2 - sin(beta2); 

d_t1 = 1/n * (t_a1 + t_b1); 
d_t2 = 1/n * (t_a2 + t_b2); 

f1 = 1 - amin / rE_mag * ( 1 - cos(dE1) ); 
f2 = 1 - amin / rE_mag * ( 1 - cos(dE2) ); 
g1 = d_t1 - 1/n * ( dE1 - sin(dE1) );  
g2 = d_t2 - 1/n * ( dE2 - sin(dE2) );  

v_E1 = 1/g1 * ( r_sunM - f1 * r_sunE );
v_E2 = 1/g2 * ( r_sunM - f2 * r_sunE );

% ellipse 1 - rv is from Sun to satellite 
rv0_e1 = [r_sunE, v_E1]; 
rv0_e2 = [r_sunE, v_E2]; 

oe_e1 = rvOrb.rv2orb(rv0_e1', mu); 
oe_e2 = rvOrb.rv2orb(rv0_e2', mu); 

%% propagate 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% Set run state 
disp('Running sim ...')

% propagate orbit 
dt = 86400/200; 
tic
[et, X_sat1] = ode45(@fn.EOM, [et_t0 : dt : et_t0 + d_t1], rv0_e1, options); 
toc 

% GMs 
mu_E = 3.986004418e5; 
mu_sun = mu; 
mu_M = 0.042828e6; 

% initialize 
a_Esat_hist = []; 
a_Msat_hist = []; 
a_sunsat_hist = []; 
a_pert_Esun_hist = []; 
a_pert_Msun_hist = []; 

const_vec = 1; 

% determine gravity 
for i = 1:length(X_sat1) 
    
    % current sun-sat vector 
    r_sunsat = X_sat1(i, 1:3);
    r_satsun = -r_sunsat; 
    
    if const_vec == 1
        r_sunE = r_sunE;   % constant vector 
        r_sunM = r_sunM; 
    else
        % Current sun-Earth vector 
        target  = 'Earth';
        X_sunE  = spice_state(et(i), target, frame, abcorr, observer); 
        r_sunE = X_sunE(1:3); 
        % Current sun-Mars vector 
        target = 'Mars'; 
        X_sunM = spice_state(et(i), target, frame, abcorr, observer); 
        r_sunM = X_sunM(1:3); 
    end 
    
        
    % Earth to satellite vector 
    r_Esun = -r_sunE;    
    r_Esat = r_Esun + r_sunsat; 

    % central body accel Earth-satellite
    a_Esat = - mu_E * r_Esat / norm(r_Esat)^3; 
    a_Esat_hist = [a_Esat_hist; norm(a_Esat)]; 
    
    % central body accel sun-satellite 
    a_sunsat = - mu_sun * r_sunsat / norm(r_sunsat)^3; 
    a_sunsat_hist = [a_sunsat_hist; norm(a_sunsat)]; 
    
    % disturbance (third body) Earth-sun 
%     a_dist = - mu_E * r_Esat / norm(r_Esat)^3 - ... 
%         mu_sun * ( r_satsun/norm(r_satsun)^3 + r_Esun/norm(r_Esun)^3); 
    % Vallado disturbance 
    a_pert_Esun = - mu_sun * ( r_satsun/norm(r_satsun)^3 - r_Esun/norm(r_Esun)^3); 
    a_pert_Esun_hist = [a_pert_Esun_hist; norm(a_pert_Esun)]; 
    
    % Mars to satellite vector 
    r_Msun = -r_sunM; 
    r_Msat = r_Msun + r_sunsat;
    
    % central body accel Mars-satellite 
    a_Msat = - mu_M * r_Msat / norm(r_Msat)^3; 
    a_Msat_hist = [a_Msat_hist; norm(a_Msat)]; 
    
    % disturbance (third body) Mars-sun 
    a_pert_Msun = - mu_sun * ( r_satsun/norm(r_satsun)^3 - r_Msun/norm(r_Msun)^3 );
    a_pert_Msun_hist = [a_pert_Msun_hist; norm(a_pert_Msun)]; 
    
end 

t_days = (et-et_t0)/86400; 
% plot accelerations 
figure
    semilogy(t_days, a_Esat_hist, 'linewidth', 2); hold on; 
    semilogy(t_days, a_sunsat_hist, '--','linewidth', 2); 
    semilogy(t_days, a_pert_Esun_hist, 'linewidth', 1.2); 
%         semilogy((et-et_t0)/86400, a_dist_Esun_P_hist, '-^'); 
    semilogy(t_days, a_Msat_hist, ':', 'linewidth', 1.2)
    semilogy(t_days, a_pert_Msun_hist, '--', 'linewidth', 1.2); 

    legend('2body E-sat', '2body sun-sat', 'pert E-sun (3rd body)', ... 
        '2body M-sat', 'pert M-sun (3rd body)', 'location', 'best')
    ylabel('km/s^2')
    xlabel('Days') 
    title('Gravitational Acceleration')
        


%% Sphere of influence 

% mass 
m_E = 5.9724e24; 
m_sun = 1988500e24; 
m_M = 0.64169e24;

% semimajor axis 
a_E = 149.598e6; 
a_M = 227.956e6; 

r_SOI_Esun = ( m_E/m_sun )^(2/5)*a_E; 
r_SOI_Msun = ( m_M/m_sun )^(2/5)*a_M; 

dt = et(2) - et(1); 

for i = 1:length(a_Msat_hist)
    
    % Earth central-perturbation body 
    ratio_Esun(i,:) = a_Esat_hist(i) / a_pert_Esun_hist(i); 
    
    % Mars central-perturbation body 
    ratio_Msun(i,:) = a_Msat_hist(i) / a_pert_Msun_hist(i); 
    
    if i > 1
        dratio_Esun(i,:) = (ratio_Esun(i,:) - ratio_Esun(i-1,:))/dt; 
        dratio_Msun(i,:) = (ratio_Msun(i,:) - ratio_Msun(i-1,:))/dt; 
    else
        dratio_Esun(i,:) = 0; 
        dratio_Msun(i,:) = 0; 
    end 
    
end 

% plot rate of ratio change 
figure
    subplot(2,1,1) 
        semilogy(t_days, dratio_Esun, 'linewidth', 2); 
    subplot(2,1,2) 
        semilogy(t_days, dratio_Msun, '--', 'linewidth', 2); 
        


%% plot 

figure()
    quiver3(0, 0, 0, r_sunE(1), r_sunE(2), r_sunE(3)); hold on; grid on; 
    quiver3(0, 0, 0, r_sunM(1), r_sunM(2), r_sunM(3))
    scatter3(0, 0, 0, 'filled')
    plot3(X_sat1(:,1), X_sat1(:,2), X_sat1(:,3))
    legend('Earth', 'Mars', 'Sun', 'Hohmann' )








