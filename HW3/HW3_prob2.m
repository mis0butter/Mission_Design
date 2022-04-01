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
T1 = T0 + dt_days1; 

disp('1st date for 60 deg phasing (JD):') 
disp(T1) 

% new time (60 deg phasing) 
[r_E1, ~, oe_E1] = xyz_ecl(T1, Earth); 
[r_M1, ~, oe_M1] = xyz_ecl(T1, Mars); 
L_E1 = oe_E1(4); 
L_M1 = oe_M1(4); 

% required time (in days) 
dL_des = 360 - 60; 
dt_days2 = (dL_des + L_E1 - L_M1) / dL_M * 100 * 365.25; 
T2 = T1 + dt_days2; 

disp('2nd date for 60 deg phasing (JD):') 
disp(T2) 

% new new time (60 deg phasing) 
[r_E2, ~, oe_E2] = xyz_ecl(T2, Earth); 
[r_M2, ~, oe_M2] = xyz_ecl(T2, Mars); 
L_E2 = oe_E2(4); 
L_M2 = oe_M2(4); 


%% now Mars to Earth - inbound conic - Part C 

% required time (in days) 
dt_days3 = (dL_des + L_M2 - L_E2) / dL_E * 100 * 365.25; 
T3 = T2 + dt_days3; 

disp('Date for 60 deg phasing (JD) inbound:') 
disp(T3) 

% new new new time (60 deg phasing) 
[r_E3, ~, oe_E3] = xyz_ecl(T3, Earth); 
[r_M3, ~, oe_M3] = xyz_ecl(T3, Mars); 
L_E3 = oe_E3(4); 
L_M3 = oe_M3(4); 


%% ok let's do a propagation ... 

T0 = 2451545.0; % units = days 

for i = 0 : 1 : dt_days1
    
    Ti = T0 + i; 
    
    % Earth 
    [r_E, r_p, oe_E] = xyz_ecl(Ti, Earth); 
    r_E_hist_T0T1(i+1,:) = r_E'; 
    oe_E_hist_T0T1(i+1,:) = oe_E'; 
    
    % Mars 
    [r_M, r_p, oe_M] = xyz_ecl(Ti, Mars); 
    r_M_hist_T0T1(i+1,:) = r_M'; 
    oe_M_hist_T0T1(i+1,:) = oe_M'; 
    
end 

for i = 0 : 1 : dt_days2
    
    Ti = T0 + dt_days1 + i; 
    
    % Earth 
    [r_E, r_p, oe_E] = xyz_ecl(Ti, Earth); 
    r_E_hist_T1T2(i+1,:) = r_E'; 
    oe_E_hist_T1T2(i+1,:) = oe_E'; 
    
    % Mars 
    [r_M, r_p, oe_M] = xyz_ecl(Ti, Mars); 
    r_M_hist_T1T2(i+1,:) = r_M'; 
    oe_M_hist_T1T2(i+1,:) = oe_M'; 
    
end 

plot3_p1p2(r_E_hist_T0T1, r_M_hist_T0T1, r_E0, r_M0, r_E1, r_M1, 'T0 to T1')

plot3_p1p2(r_E_hist_T1T2, r_M_hist_T1T2, r_E1, r_M1, r_E2, r_M2, 'T1 to T2')


%% some tests 

r_dep = r_E1; 
r_arr = r_M2; 

r_dep = r_dep / norm(r_dep); 
r_arr = r_arr / norm(r_arr); 

acosd(dot(r_arr, r_dep))

%% subfunctions 

function h = plot3_quiver(r1, r2, style)

    h = quiver3(r1(1), r1(2), r1(3), r2(1), r2(2), r2(3), style); 

end 

function plot3_p1p2(r_E_hist, r_M_hist, r_E0, r_M0, r_E1, r_M1, ftitle)

figure('name', ftitle)

    fn.plot3_xyz(r_E_hist, 'b--'); 
    hold on; grid on; 
    fn.plot3_xyz(r_E_hist(1:100,:), 'b', 2); 
    fn.plot3_xyz(r_E_hist(100,:), 'b^'); 
    fn.plot3_xyz(r_M_hist, 'r--'); 
    fn.plot3_xyz(r_M_hist(1:100,:), 'r', 2); 
    fn.plot3_xyz(r_M_hist(100,:), 'r^'); 

    % E0 point 
    fn.plot3_xyz(r_E0', 'bo'); 
    plot3_quiver([0 0 0], r_E0, 'b'); 
    txt = '    E_i'; 
    text(r_E0(1), r_E0(2), r_E0(3), txt)

    % M0 point 
    fn.plot3_xyz(r_M0', 'ro'); 
    plot3_quiver([0 0 0], r_M0, 'r'); 
    txt = '    M_i'; 
    text(r_M0(1), r_M0(2), r_M0(3), txt)

    % E1 point 
    fn.plot3_xyz(r_E1', 'bp'); 
    plot3_quiver([0 0 0], r_E1, 'b'); 
    txt = '    E_f'; 
    text(r_E1(1), r_E1(2), r_E1(3), txt)

    % M1 point 
    fn.plot3_xyz(r_M1', 'rp'); 
    plot3_quiver([0 0 0], r_M1, 'r'); 
    txt = '    M_f'; 
    text(r_M1(1), r_M1(2), r_M1(3), txt)
    
    xlim([-2 2])
    ylim([-2 2])
%     axis equal 
    xlabel('x (AU)') 
    ylabel('y (AU)') 
    zlabel('z (AU)') 
    view(0,90)

    sgtitle(ftitle)

end 



