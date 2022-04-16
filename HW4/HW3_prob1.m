%% HW 2 
% Junette Hsin 

% Keplerian elements 
% a     = semi-major axis 
% e     = eccentricity 
% i     = inclination 
% L     = mean longitude 
% wbar  = longitude of perihelion 
% Omega = longitude of ascending node 

clear; clc 

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
Mars.L0 =     -4.55343205 + 360; 
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

%% find launch date 

T0 = 2451545.0; % units = days 

T = T0;
[dt_days1, tof1_1, tof1_2] = launch_date(T, Earth, Mars, mu); 
% find 1st launch date 
while abs(dt_days1 - tof1_1) > 1 && abs(dt_days1 - tof1_2) > 1

    T = T + 1; 
    [dt_days1, tof1_1, tof1_2] = launch_date(T, Earth, Mars, mu); 
    
end 

% check 
[r_E0, ~, oe_E0] = xyz_ecl(T, Earth); 
[r_M0, ~, oe_M0] = xyz_ecl(T, Mars); 
[r_E1, ~, oe_E1] = xyz_ecl(T + dt_days1, Earth); 
[r_M1, ~, oe_M1] = xyz_ecl(T + dt_days1, Mars); 
phi = acosd(dot(r_E0, r_M1) / ( norm(r_E0)*norm(r_M1) ))
disp('1st launch window: ')
disp('E0 longitude = ') 
L_E0 = oe_E0(4)
disp('M0 longitude = ') 
L_M1 = oe_M0(4)
disp('M1 longitude = ') 
L_M1 = oe_M1(4)

T0 = T; 
T1 = T0 + dt_days1; 

% ok let's do a propagation ... 

r_E_hist_T0T1 = []; 
r_M_hist_T0T1 = []; 
for i = 0 : 0.1 : dt_days1
    
    Ti = T0 + i; 
    
    % Earth 
    [r_E, r_p, oe_E] = xyz_ecl(Ti, Earth); 
    r_E_hist_T0T1 = [r_E_hist_T0T1; r_E']; 
%     oe_E_hist_T0T1(i+1,:) = oe_E'; 
    
    % Mars 
    [r_M, r_p, oe_M] = xyz_ecl(Ti, Mars); 
    r_M_hist_T0T1 = [r_M_hist_T0T1; r_M']; 
%     oe_M_hist_T0T1(i+1,:) = oe_M'; 
    
end 

txt = {sprintf('T0 = %.10g JD', T0); ... 
    sprintf('T1 = %.10g JD', T1); ... 
    sprintf('TOF = %.10g days', dt_days1) ... 
    }; 
plot3_p1p2(r_E_hist_T0T1, r_M_hist_T0T1, r_E0, r_M0, r_E1, r_M1, 'Prob 1 - 1st Launch Window', txt)
fn.save_pdf(gcf)


%% find 2nd launch date 

T = T1;
[dt_days1, tof1_1, tof1_2] = launch_date(T, Earth, Mars, mu); 
% find 1st launch date 
while abs(dt_days1 - tof1_1) > 1 && abs(dt_days1 - tof1_2) > 1

    T = T + 1; 
    [dt_days1, tof1_1, tof1_2] = launch_date(T, Earth, Mars, mu); 
    
end 

% check 
[r_E0, ~, oe_E0] = xyz_ecl(T, Earth); 
[r_M0, ~, oe_M0] = xyz_ecl(T, Mars); 
[r_E1, ~, oe_E1] = xyz_ecl(T + dt_days1, Earth); 
[r_M1, ~, oe_M1] = xyz_ecl(T + dt_days1, Mars); 
phi = acosd(dot(r_E0, r_M1) / ( norm(r_E0)*norm(r_M1) ))
disp('2nd launch window: ')
disp('E0 longitude = ') 
L_E0 = oe_E0(4)
disp('M0 longitude = ') 
L_M1 = oe_M0(4)
disp('M1 longitude = ') 
L_M1 = oe_M1(4)

T0 = T; 
T1 = T0 + dt_days1; 

% ok let's do a propagation ... 

r_E_hist_T0T1 = []; 
r_M_hist_T0T1 = []; 
for i = 0 : 0.1 : dt_days1
    
    Ti = T0 + i; 
    
    % Earth 
    [r_E, r_p, oe_E] = xyz_ecl(Ti, Earth); 
    r_E_hist_T0T1 = [r_E_hist_T0T1; r_E']; 
%     oe_E_hist_T0T1(i+1,:) = oe_E'; 
    
    % Mars 
    [r_M, r_p, oe_M] = xyz_ecl(Ti, Mars); 
    r_M_hist_T0T1 = [r_M_hist_T0T1; r_M']; 
%     oe_M_hist_T0T1(i+1,:) = oe_M'; 
    
end 

txt = {sprintf('T0 = %.10g JD', T0); ... 
    sprintf('T1 = %.10g JD', T1); ... 
    sprintf('TOF = %.10g days', dt_days1) ... 
    }; 
plot3_p1p2(r_E_hist_T0T1, r_M_hist_T0T1, r_E0, r_M0, r_E1, r_M1, 'Prob 1 - 2nd Launch Window', txt)
fn.save_pdf(gcf)


%% inbound conic 

T = T1;
[dt_days1, tof1_1, tof1_2] = launch_date(T, Mars, Earth, mu); 
% find 1st launch date 
while abs(dt_days1 - tof1_1) > 1 && abs(dt_days1 - tof1_2) > 1

    T = T + 1; 
    [dt_days1, tof1_1, tof1_2] = launch_date(T, Mars, Earth, mu); 
    
end 

% check 
[r_E0, ~, oe_E0] = xyz_ecl(T, Earth); 
[r_M0, ~, oe_M0] = xyz_ecl(T, Mars); 
[r_E1, ~, oe_E1] = xyz_ecl(T + dt_days1, Earth); 
[r_M1, ~, oe_M1] = xyz_ecl(T + dt_days1, Mars); 
phi = acosd(dot(r_E0, r_M1) / ( norm(r_E0)*norm(r_M1) ))
disp('Inbound conic launch window: ')
disp('E0 longitude = ') 
L_E0 = oe_E0(4)
disp('E1 longitude = ') 
L_E1 = oe_E1(4)
disp('M0 longitude = ') 
L_M1 = oe_M0(4)

T0 = T; 
T1 = T0 + dt_days1; 

% ok let's do a propagation ... 

r_E_hist_T0T1 = []; 
r_M_hist_T0T1 = []; 
for i = 0 : 0.1 : dt_days1
    
    Ti = T0 + i; 
    
    % Earth 
    [r_E, r_p, oe_E] = xyz_ecl(Ti, Earth); 
    r_E_hist_T0T1 = [r_E_hist_T0T1; r_E']; 
%     oe_E_hist_T0T1(i+1,:) = oe_E'; 
    
    % Mars 
    [r_M, r_p, oe_M] = xyz_ecl(Ti, Mars); 
    r_M_hist_T0T1 = [r_M_hist_T0T1; r_M']; 
%     oe_M_hist_T0T1(i+1,:) = oe_M'; 
    
end 

txt = {sprintf('T0 = %.10g JD', T0); ... 
    sprintf('T1 = %.10g JD', T1); ... 
    sprintf('TOF = %.10g days', dt_days1) ... 
    }; 
plot3_p1p2(r_E_hist_T0T1, r_M_hist_T0T1, r_E0, r_M0, r_E1, r_M1, 'Prob 1 - Inbound Conic Launch Window', txt)
fn.save_pdf(gcf)


%% subfunctions 

function h = plot3_quiver(r1, r2, style)

    h = quiver3(r1(1), r1(2), r1(3), r2(1), r2(2), r2(3), style); 

end 

function plot3_p1p2(r_E_hist, r_M_hist, r_E0, r_M0, r_E1, r_M1, ftitle, plt_txt)

if ~exist('plt_txt', 'var')
    plt_txt = ''; 
end 

figure('name', ftitle, 'position', [100 100 600 700])

    subplot(3,1,1:2)

        fn.plot3_xyz(r_E_hist, 'b--'); 
        hold on; grid on; 
        n = 1000; 
        fn.plot3_xyz(r_E_hist(1:n,:), 'b', 2); 
        fn.plot3_xyz(r_E_hist(n,:), 'b^'); 
        fn.plot3_xyz(r_M_hist, 'r--'); 
        fn.plot3_xyz(r_M_hist(1:n,:), 'r', 2); 
        fn.plot3_xyz(r_M_hist(n,:), 'r^'); 

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
        
    subplot(3,1,3) 
    
        fn.plt_txt(plt_txt) 

    sgtitle(ftitle)

end 

function [dt_days1, tof1, tof2] = launch_date(T0, dep, arr, mu)


    [r_dep0, ~, oe_dep0] = xyz_ecl(T0, dep); 
    [r_arr0, ~, oe_arr0] = xyz_ecl(T0, arr); 

    % km to au 
    km2au = 6.6845871226706E-9; 
    au2km = 1/km2au; 

    % mean longitude for Departure and Arrival
    L_dep0 = oe_dep0(4); 
    L_arr0 = oe_arr0(4); 

    % desired delta longitude 
    dL_des = 60; 

    % required time (in days) 

%     dL = dL_des + L_M0 - L_E0; 
%     if dL < 0
%         dL = dL + 360; 
%     elseif dL > 360 
%         while dL > 360 
%             dL = dL - 360; 
%         end 
%     end 
    
    % desired longitude for arrival 
    L_des = L_dep0 + dL_des; 
    
    % delta longitude 
    dL = L_des - L_arr0; 
    if dL < 0 
        dL = dL + 360; 
    elseif dL > 360 
        dL = dL - 360; 
    end
    
    dL_arr = arr.dL;
    dt_days1 = dL / dL_arr * 100 * 365.25; 
%     dt_days1 = (dL_des + L_E0 - L_M0) / dL_M * 100 * 365.25; 
    T1 = T0 + dt_days1; 

    % new time (60 deg phasing) 
    [r_dep1, ~, oe_dep1] = xyz_ecl(T1, dep); 
    [r_arr1, ~, oe_arr1] = xyz_ecl(T1, arr); 
    L_dep1 = oe_dep1(4); 
    L_arr1 = oe_arr1(4); 

    [ell_1, ell_2] = bett_lambert(r_dep0' * au2km, r_arr1' * au2km, mu); 
    tof1 = ell_1.tof / 86400; 
    tof2 = ell_2.tof / 86400; 

end 



