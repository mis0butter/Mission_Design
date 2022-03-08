% 3. [Problem 6.3, Capderou] Calculate the dates during the year 1999 for which the local mean
% time of the ascending node crossing is the same for the satellites TRMM and Resurs-O1-4.
% TRMM flew in a near-circular 350 km altitude orbit at 35° inclination. TRMM crossed the
% ascending node at time 1999-01-21 20:43:47 (UT) at geographic longitude of 5.157° West.
% The satellite Resurs-O1-4 flew in a sun-synchronous orbit at 22:20 local mean time. Show all
% calculations.

% [NOTE: To quote from Capderou: “In order to study the Earth’s radiation budget, TRMM
% and Resurs-O1-4 were equipped with the CERES and ScaRaB instruments, respectively.
% A joint measurement campaign was organized in January and February 1999. The aim
% was to compare the measurements obtained for the same region viewed by the two
% instruments at roughly the same time (with a leeway of ± 15 minutes).”
% The solution to this problem has been diagrammed in a previous email – but it is
% important to work through the numbers, nevertheless.]

%% 3. [Problem 6.3, Capderou] 

% Calculate the dates during the year 1999 for which the local mean time 
% of the ascending node crossing is the same for the satellites TRMM and 
% Resurs-O1-4.

% TRMM flew in a near-circular 350 km altitude orbit at 35° inclination. 
% TRMM crossed the ascending node at time 1999-01-21 20:43:47 (UT) at 
% geographic longitude of 5.157° West. latitude = 0 deg 

global mu  
global eop_data 
eop_data = load('finals_iau1980.txt'); 

% Initial conditions 
mu = 3.986004415e14; %m^3/s^2

% a = Earth radius + 350 km   
%   Earth radius = 6378136.3 m 
a = 6.7281363e6;  % m 
% i = 35 deg 
i = 35 * pi/180; 
e = 0.0001; 
% e = 0.000; 
% LAN = 5.157 deg W longitude @ 1999-01-21 20:43:47 (UT)
LAN = 5.157 * pi/180; 
% w = 0
w = 0; 
% true anomaly (angle between periapsis and current orbit position) = 

oe0_TRM = [a; e; i; w; LAN; LAN]; 

rv0_TRM = oe2rv(0, oe0_TRM); 


%% TRMM flew in a near-circular 350 km altitude orbit at 35° inclination. 
% TRMM crossed the ascending node at time 1999-01-21 20:43:47 (UT) at 
% geographic longitude of 5.157° West. latitude = 0 deg 

% LAT = 0, LON = -5.157, ALT = 350e3 M 
% X : 6700.902   km
% Y : -604.76   km
% Z : 0   km

% JD time (1999-01-21 20:43:47 UTC)
JD = 2451200.36389; 

% ECEF frame 
r_ECEF = [ 6700.902; -604.76; 0 ] * 1000; % m 
[r_ECI] = fn.ECEFtoECI_r(JD, r_ECEF)

% Keplerian orbit  
v_mag = sqrt(mu/norm(r_ECEF)); 

% DCM between ECEF and 35 deg tilt 
ECEF_DCM_tilt35 = cosd([0, 90, 90; 90, 35, 125; 90, 55, 35]); 

% z vector of orbit plane 
orbit_n = ECEF_DCM_tilt35(:, end); 

% get velocity vector 
v_vec = cross(orbit_n, r_ECEF / norm(r_ECEF)); 
v_ECEF = v_vec * v_mag; 
v_ECI = fn.ECEFtoECI_r(JD, v_ECEF); 

% initial position and velocity 
rv0 = [r_ECI; v_ECI]; 

% PROPAGATE ORBIT 
t0 = 0; 
dt = 200; 
tf = 60 * 60 * 24 * 365; 

% Options and parameters 
toler   = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
options = odeset('reltol', toler, 'abstol', toler ); 

% Solve ODE 
[t,x] = ode45(@TwoBod_6states, [t0 : dt : tf], [rv0], options); 
r_x  = x(:, 1); r_y  = x(:, 2); r_z  = x(:, 3); 

% Plot 
figure('name', 'Problem 0: Relative 2 Body EOM') 
    plot3(r_x, r_y, r_z, '-o', 'markersize', 2); 
    hold on; grid on; 
    scatter3(r_x(1), r_y(1), r_z(1)); 
    scatter3(r_x(end), r_y(end), r_z(end), 'kx'); 
    plot3(0, 0, 0, 'o', 'markersize', 6, 'linewidth', 2)
    
    legend('orbit', 'start', 'finish', 'Earth')
    title('Problem 0: Relative 2 Body EOM') 
    xlabel('x (LU)'); ylabel('y (LU)'); zlabel('z (LU)'); 
    
%% convert to orbital elements 



for i = 1:length(x)
    oe_TRM(i,:) = rv2oe(x(i, 1:6)); 
end 


%% 

% The satellite Resurs-O1-4 flew in a sun-synchronous orbit at 22:20 local 
% mean time. Show all calculations. Almost circular? 
% 1999-01-21 22:20 local mean time --> Jan 22, 1999 4:20 UTC 

% h = 820 km (altitude) 
%   a = 6.7281363e6 + 820000 m 
a = 7.5481363e+6; 
% i = 98.8 deg 
i = 98.8 * pi/180; 
e = 0.0001; 
% LAN = 5.157 deg W longitude @ 1999-01-21 20:43:47 (UT)
LAN = 5.157 * pi/180; 
% w = 0; 
w = 0; 

oe0_RSS = [a; i; e; w; LAN; LAN]; 

rv0_RSS = oe2rv(0, oe0_RSS); 


%% Propagate orbit 

t0 = 0; 
dt = 100; 
tf = 60 * 60 * 24 * 365; 

% Options and parameters 
toler   = 1e-8;         % 1e-14 accurate; 1e-6 coarse 
options = odeset('reltol', toler, 'abstol', toler ); 

% Solve ODE 
[t,x] = ode45(@TwoBod_6states, [t0 : dt : tf], [rv0_TRM], options); 
r_x  = x(:, 1); r_y  = x(:, 2); r_z  = x(:, 3); 

% Plot 
figure('name', 'Problem 0: Relative 2 Body EOM') 
    plot3(r_x, r_y, r_z, '-o', 'markersize', 2); 
    hold on; grid on; 
    scatter3(r_x(1), r_y(1), r_z(1)); 
    scatter3(r_x(end), r_y(end), r_z(end), 'kx'); 
    plot3(0, 0, 0, 'o', 'markersize', 6, 'linewidth', 2)
    
    legend('orbit', 'start', 'finish', 'Earth')
    title('Problem 0: Relative 2 Body EOM') 
    xlabel('x (LU)'); ylabel('y (LU)'); zlabel('z (LU)'); 


%% Convert between Earth longitude and J2000 

addpath(genpath('mice')) 
addpath(genpath('spice_data'))

clear; clc;

%  Load kernel file 
cspice_furnsh( 'spice_data/naif0011.tls' )
cspice_furnsh( 'spice_data/de421.bsp' )       
cspice_furnsh( 'spice_data/pck00010.tpc ') 

%  Define parameters for a state lookup:
% t0      = 'Oct 20, 2020 11:00 AM CST'; 
t0      = 'Jan 21, 1999 20:44 UTC'; 
abcorr  = 'NONE';

%  Convert the epoch to ephemeris time. 
et_t0   = cspice_str2et( t0 );

% Earth wrt Moon 
target      = 'Earth';
frame       = 'J2000';
observer    = 'Moon';
days        = 0; 

% CSPICE_PXFORM returns the matrix that transforms position
% vectors from one specified frame to another at a specified epoch.
% Retrieve the transformation matrix FROM (1) TO (2) at epoch (3).
%         mat{i,:}    = cspice_pxform( 'IAU_MERCURY', 'ECLIPJ2000', epoch );
mat    = cspice_pxform( 'IAU_EARTH', 'IAU_MOON', et_t0);

%  Look-up the state for the defined parameters.
starg  = mice_spkezr( target, et_t0, frame, abcorr, observer);
rv     = starg.state(1:6); 



%% 







