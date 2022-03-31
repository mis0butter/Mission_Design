%% HW 2 
% Junette Hsin 

% Keplerian elements 
% a     = semi-major axis 
% e     = eccentricity 
% i     = inclination 
% L     = mean longitude 
% wbar  = longitude of perihelion 
% Omega = longitude of ascending node 

close all; clear; clc 

% Earth 
Earth.a0 =      1.00000261; 
Earth.da =      0.00000562; 
Earth.e0 =      0.01671123; 
Earth.de =     -0.00004392; 
Earth.I0 =     -0.00001531; 
Earth.dI =     -0.01294668;
Earth.L0 =    100.46457166; 
Earth.dL =  35999.37244981; 
Earth.wbar0 = 102.93768193; 
Earth.dwbar =   0.32327364;
Earth.Omega0 =  0; 
Earth.dOmega =  0; 

% Mars 
Mars.a0 =      1.52371034; 
Mars.da =      0.00001847; 
Mars.e0 =      0.09339410; 
Mars.de =      0.00007882; 
Mars.I0 =      1.84969142; 
Mars.dI =     -0.00813131; 
Mars.L0 =     -4.55343205; 
Mars.dL =  19140.30268499; 
Mars.wbar0 = -23.94362959; 
Mars.dwbar =   0.44441088; 
Mars.Omega0 = 49.55953891; 
Mars.dOmega = -0.29257343; 

% sun mu (m^3/s^2)
mu_sun_m = 1.32712440018e20; 
mu_sun_km = mu_sun_m / (1000^3); 
mu = mu_sun_km; 

%% synodic period 

dL_E = Earth.dL; 
dL_M = Mars.dL;         % degrees per century 

% period of Earth 
T_E = 365.25; 

% period of Mars 
dL_M_degyear = dL_M / 100;      % degs per year
T_M_years = (1/dL_M_degyear) * 360;            % period (days per 360 deg) 

% synodic period 
SP_M = T_M_years / abs(T_M_years - 1); 

disp('Synodic period:') 
disp(SP_M) 

%% phase angle 60 deg - Part A and B 

T0 = 2451545.0; % units = days 
[r_E0, ~, oe_E0] = xyz_ecl(T0, Earth); 
[r_M0, ~, oe_M0] = xyz_ecl(T0, Mars); 

% km to au 
km2au = 6.6845871226706E-9; 
au2km = 1/km2au; 

% mean longitude for Earth and Mars 
L_E0 = oe_E0(4); 
L_M0 = oe_M0(4); 

% desired delta longitude 
dL_des = 60; 

% required time (in days) 
dt_days1 = (dL_des + L_E0 - L_M0) / dL_M * 100 * 365.25; 

% new time (60 deg phasing) 
T1 = T0 + dt_days; 
[r_E1, ~, oe_E1] = xyz_ecl(T1, Earth); 
[r_M1, ~, oe_M1] = xyz_ecl(T1, Mars); 
L_E1 = oe_E1(4); 
L_M1 = oe_M1(4); 

% required time (in days) 
dt_days2 = (dL_des + L_E1 - L_M1) / dL_M * 100 * 365.25; 

% new new time (60 deg phasing) 
T2 = T1 + dt_days2; 
[r_E2, ~, oe_E2] = xyz_ecl(T2, Earth); 
[r_M2, ~, oe_M2] = xyz_ecl(T2, Mars); 
L_E2 = oe_E2(4); 
L_M2 = oe_M2(4); 


%% now Mars to Earth - inbound conic - Part C 

% required time (in days) 
dt_days3 = (dL_des + L_M2 - L_E2) / dL_E * 100 * 365.25; 

% new new new time (60 deg phasing) 
T3 = T2 + dt_days3; 
[r_E3, ~, oe_E3] = xyz_ecl(T3, Earth); 
[r_M3, ~, oe_M3] = xyz_ecl(T3, Mars); 
L_E3 = oe_E3(4); 
L_M3 = oe_M3(4); 






