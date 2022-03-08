%% HW 2 
% Junette Hsin 

% Keplerian elements 
% a     = semi-major axis 
% e     = eccentricity 
% i     = inclination 
% L     = mean longitude 
% wbar  = longitude of perihelion 
% Omega = longitude of ascending node 

% close all; 
clear; 

%% transfer angle = 75 deg 

addpath(genpath('mice')); 
addpath(genpath('spice_data')); 

%  Load kernel file 
cspice_furnsh( 'spice_data/naif0011.tls' )
cspice_furnsh( 'spice_data/de421.bsp' )       
cspice_furnsh( 'spice_data/pck00010.tpc ') 

%  Define parameters for a state lookup:
% t0      = 'Oct 20, 2020 11:00 AM CST'; 
t0      = 'May 22, 2000'; 
abcorr  = 'NONE';

%  Convert the epoch to ephemeris time (secs) 
et_t0   = cspice_str2et( t0 );

% get states --> Sun to Earth 
target   = 'Earth';
frame    = 'J2000';
observer = 'Sun';
abcorr   = 'NONE';

% get sun position 
et = et_t0;    % propagate ephemeris time by 1 day in secs 
X_sunE  = spice_state(et, target, frame, abcorr, observer); 

% get states --> Sun to Mars
target   = 'Mars';
frame    = 'J2000';
observer = 'Sun';
abcorr   = 'NONE';

% get sun position 
et = et_t0;    % propagate ephemeris time by 1 day in secs 
X_sunM  = spice_state(et, target, frame, abcorr, observer); 

% get angle between Earth and Mars velocities 
r_E = X_sunE(1:3); 
r_M = X_sunM(1:3); 
v_E = X_sunE(4:6); 
v_M = X_sunM(4:6); 

% angle and velocity angles 
phi_r = acosd( dot(r_E, r_M) / (norm(r_E)*norm(r_M)) ); 
phi_v = acosd( dot(v_E, v_M) / (norm(v_E)*norm(v_M)) ); 

%% phi = 75 

phi_des = 75; 
lambert_prob 

%% phi = 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180 degrees

for phi = 30 : 30 : 180 
    
    phi_des = phi; 
    lambert_prob 
    
end 












