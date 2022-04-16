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
target   = 'Sun';
frame    = 'J2000'; 
observer = 'Earth';
abcorr   = 'NONE'; 
et = et_t0;    

% get Earth-sun position 
X_Esun  = spice_state(et, target, frame, abcorr, observer); 

% get Earth-Mars position 
target = 'Mars'; 
X_EM = spice_state(et, target, frame, abcorr, observer); 

% get angle between Earth and Mars velocities 
r_Esun = X_Esun(1:3); 
r_EM   = X_EM(1:3); 
v_Esun = X_Esun(4:6); 
v_EM   = X_EM(4:6); 

% position and velocity angles 
phi_r = acosd( dot( r_Esun, r_EM ) / ( norm(r_Esun)*norm(r_EM)) ); 
phi_v = acosd( dot( v_Esun, v_EM ) / ( norm(v_Esun)*norm(v_EM)) ); 

% desired phase angle 
phi_des = 0.7073; 

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
    X_EM = spice_state(et, target, frame, abcorr, observer); 
    r_EM = X_EM(1:3); 
    v_EM = X_EM(4:6); 
    
    % get angle 
    phi_r = acosd( dot( r_Esun, r_EM ) / ( norm(r_Esun)*norm(r_EM)) )
    phi_v = acosd( dot( v_Esun, v_EM ) / ( norm(v_Esun)*norm(v_EM)) ); 
    
end 

%% bettadpur lambert 

% Lambert - from VALLADO ed. 4, pgs 467 - 475 

% Arrival (Mars), AU units 
rM_mag = norm(r_EM); 
rsun_mag = norm(r_Esun); % Vallado method ... 
cos_dv = dot(r_Esun, r_EM) / (rsun_mag * rM_mag); 

% chord 
c = sqrt( rsun_mag^2 + rM_mag^2 - 2*rsun_mag*rM_mag*cos_dv ); 

% semiperimeter 
s = ( rM_mag + rsun_mag + c ) / 2; 

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

f1 = 1 - amin / rsun_mag * ( 1 - cos(dE1) ); 
f2 = 1 - amin / rsun_mag * ( 1 - cos(dE2) ); 
g1 = d_t1 - 1/n * ( dE1 - sin(dE1) );  
g2 = d_t2 - 1/n * ( dE2 - sin(dE2) );  

v_E1 = 1/g1 * ( r_EM - f1 * r_Esun );
v_E2 = 1/g2 * ( r_EM - f2 * r_Esun ); 

% ellipse 1 - rv is from Earth to satellite 
rv0_e1 = [r_Esun, v_E1]; 
oe_e1 = rvOrb.rv2orb(rv0_e1', mu); 

% propagate 

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





%% plot 

figure()
    quiver3(0, 0, 0, r_Esun(1), r_Esun(2), r_Esun(3)); hold on; grid on; 
    quiver3(0, 0, 0, r_EM(1), r_EM(2), r_EM(3))
    scatter3(0, 0, 0, 'filled')
    plot3(X_sat1(:,1), X_sat1(:,2), X_sat1(:,3))
    legend('Sun', 'Mars', 'Earth', 'Hohmann' )


















