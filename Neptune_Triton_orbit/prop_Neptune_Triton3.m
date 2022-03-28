% HW 5 
% Junette Hsin 

% clear; clc 
addpath(genpath('mice')); 
addpath(genpath('spice_data')); 

% Load SPICE kernel file 
cspice_furnsh( 'spice_data/naif0011.tls' )
cspice_furnsh( 'spice_data/de421.bsp' )       
cspice_furnsh( 'spice_data/pck00010.tpc ') 
cspice_furnsh( 'spice_data/nep095.bsp' )

format long g 

% constants 
constants 


%% propagate Triton 

% initial state initial guess (M --> KM)

t0 = 'May 1, 2049, 00:00:00 UTC'; 

% integrate Triton 
abcorr  = 'NONE';

%  Convert the epoch to ephemeris time (secs) 
et_t0   = cspice_str2et( t0 );

% get states --> Neptune to Triton
target   = 'Triton';
frame    = 'ECLIPJ2000';
observer = 'Neptune';
abcorr   = 'NONE';

X_NT = []; 
OE_T = []; 

T0_T = 507772.8; 
dt = 10; 
for i = 0 : dt : T0_T - 86400
    
    % propagate by 0.1 day 
    et = et_t0 + i; 
    
    % get velocity 
    X  = spice_state(et, target, frame, abcorr, observer); 
    r_T = X(1:3); 
    v_T = X(4:6); 
    
    OE = rvOrb.rv2orb(X, const.muN); 
    
    X_NT = [X_NT; X]; 
    OE_T = [OE_T; OE]; 
    
end 

%% initial state 

% hyperbolic orbit 
% rv0_test = [19813.3, -16908.2, 2612.7, -22.953, -13.324, 11.316]; 

% "beginning" of hyperbolic orbit 
rv0_sat = [ 
          610470.376325053
          739754.071895096
         -468914.276897989
         -9.72778840140298
         -12.8560823075242
          7.93760951868089]; 

OE0_sat = rvOrb.rv2orb(rv0_sat, const.muN); 

% satellite period 
T_sat = 2*pi*sqrt( OE0_sat(1)^3 / const.muN ); 


%% integrate EOM - satellite 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% Set run state 
disp('Running sim ...')

dt = 10; 

% propagate hyperbolic orbit 
tic
[t, X_Nsat_h] = ode45(@fn.EOM, [0 : dt : T0_T/4], rv0_sat, options); 
toc 
disp('Pos and Vel end: ')
disp(X_Nsat_h(end, 1:6)'); 

% ------------------------------------------------------------------------ 
% 1st aerocapture ellipse (decrease apoapsis) 

% compute rnorm 
rnorm_Nsat = []; 
for i = 1:length(X_Nsat_h)
    rnorm_Nsat(i,:) = norm(X_Nsat_h(i,1:3)); 
end 

% find index of min norm - save state at periapsis 
i_min = find(rnorm_Nsat == min(rnorm_Nsat)); 
X_Nsat_min = X_Nsat_h(i_min, :); 

% augment OE "delta v" aerocapture from hyperbolic --> elliptical 
OE_Nsat_min = rvOrb.rv2orb(X_Nsat_min, const.muN); 
OE_Nsat_min(1) = OE0_T(1); 
OE_Nsat_min(2) = 0.9; 

% obtain state of augmented elliptical orbit 
rv0_test = rvOrb.orb2rv(OE_Nsat_min, const.muN); 

% extract velocity direction, save over periapsis state 
X_Nsat_min(4:6) = rv0_test(4:6)*1.2; 
OE_Nsat_min = rvOrb.rv2orb(X_Nsat_min, const.muN); 
T_Nsat_min = 2*pi*sqrt( OE_Nsat_min(1)^3 / const.muN ); 

% propagate elliptical orbit 
[t, X_Nsat_e1] = ode45(@fn.EOM, [0: dt : T_Nsat_min], X_Nsat_min, options); 

% ------------------------------------------------------------------------ 
% 2nd aerocapture ellipse (decrease apoapsis) 

% obtain state of augmented elliptical orbit 
rv0_test = rvOrb.orb2rv(OE_Nsat_min, const.muN); 

% extract velocity direction, save over periapsis state 
X_Nsat_min(4:6) = rv0_test(4:6)*0.98; 
OE_Nsat_min = rvOrb.rv2orb(X_Nsat_min, const.muN); 
T_Nsat_min = 2*pi*sqrt( OE_Nsat_min(1)^3 / const.muN ); 

% propagate elliptical orbit 
[t, X_Nsat_e2] = ode45(@fn.EOM, [0: dt : T_Nsat_min*0.9], X_Nsat_min, options); 

% ------------------------------------------------------------------------ 
% raise periapsis  

% compute rnorm again 
rnorm_Nsat = []; 
for i = 1:length(X_Nsat_e2)
    rnorm_Nsat(i,:) = norm(X_Nsat_e2(i,1:3)); 
end 

% now find index of MAX norm - save state at apoapsis 
i_max = find(rnorm_Nsat == max(rnorm_Nsat)); 
X_Nsat_max = X_Nsat_e2(i_max, :); 

% augment OE "delta v" raise periapsis 
OE_Nsat_max = rvOrb.rv2orb(X_Nsat_max, const.muN); 
OE_Nsat_max(1) = OE_Nsat_max(1)*0.6; 

% obtain state of raised periapsis elliptical orbit 
rv0_test = rvOrb.orb2rv(OE_Nsat_max, const.muN); 
X_Nsat_max(4:6) = rv0_test(4:6); 
OE_Nsat_max = rvOrb.rv2orb(X_Nsat_max, const.muN); 
T_Nsat_max = 2*pi*sqrt( OE_Nsat_max(1)^3 / const.muN ); 

% propagate raised periapsis elliptical orbit 
[t, X_Nsat_ep] = ode45(@fn.EOM, [0: dt : T_Nsat_max * 0.8], X_Nsat_max, options); 

% ------------------------------------------------------------------------ 
% 1st inclination change maneuver 

[t, X_Nsat_ei1] = incl_change(X_Nsat_ep, OE_T, dt, const, options); 
[t, X_Nsat_ei2] = incl_change(X_Nsat_ei1, OE_T, dt, const, options); 
[t, X_Nsat_ei3] = incl_change(X_Nsat_ei2, OE_T, dt, const, options); 
[t, X_Nsat_ei4] = incl_change(X_Nsat_ei3, OE_T, dt, const, options); 


%% plot 

% Neptune sphere 
[X, Y, Z] = sphere; 
XN = X * const.RN; YN = Y * const.RN; ZN = Z * const.RN; 

% vernal equinox 
v_eq = [rnorm_T, 0, 0]; 


plot_orbit = 1; 
if plot_orbit == 1
    ftitle = 'Orbit around Neptune'; 
    figure('name', ftitle); 
    
        % hyperbolic, 1st aerocapture, 2nd aerocapture 
        subplot(2,2,1) 
            plot_NT(const, X_NT, rnorm_T)

            % hyperbolic trajectory 
            c = [1 0 1]; 
            plot_traj(X_Nsat_h, c); 
            plot3(X_Nsat_min(1), X_Nsat_min(2), X_Nsat_min(3), 'mp'); 

            % 1st elliptical aerocapture 
            c = [0 1 1]; 
            plot_traj(X_Nsat_e1, c); 

            % 2nd elliptical aerocapture 
            plot_traj(X_Nsat_e2, c); 

            % APOAPSIS 
            plot3(X_Nsat_e2(i_max,1), X_Nsat_e2(i_max,2), X_Nsat_e2(i_max,3), 'p'); 
            
            title('Hyperb, 1st aero, 2nd aero') 
            
        % 2nd aerocapture, raise periapsis 
        subplot(2,2,2) 
            plot_NT(const, X_NT, rnorm_T)

            % 2nd elliptical aerocapture 
            c = [0 1 1]; 
            plot_traj(X_Nsat_e2, c); 
            
            % APOAPSIS 
            plot3(X_Nsat_e2(i_max,1), X_Nsat_e2(i_max,2), X_Nsat_e2(i_max,3), 'p'); 

            % raise periapsis 
            c = [0.8500 0.3250 0.0980]; 
            plot_traj(X_Nsat_ep, c)
            
            title('2nd aero, raise peri') 
            
        % raise periapsis, inclination change 
        subplot(2,2,3) 
            plot_NT(const, X_NT, rnorm_T)

            % raise periapsis 
            c = [0.8500 0.3250 0.0980]; 
            plot_traj(X_Nsat_ep, c)

            % APOAPSIS 
            plot3(X_Nsat_e2(i_max,1), X_Nsat_e2(i_max,2), X_Nsat_e2(i_max,3), 'p'); 
        
            % RAAN line node 
            a = OE_Nsat_min(1); 
            quiver3(0, 0, 0, n(1)*a, n(2)*a, n(3)*a, 'k', 'linewidth', 1.5); 
            plot3(rv_di(1), rv_di(2), rv_di(3), 'kp')
                txt = 'line node'; 
                text(n(1)*a, n(2)*a, n(3)*a, txt)
            quiver3(0, 0, 0, -n(1)*a*2, -n(2)*a*2, -n(3)*a*2, 'k', 'linewidth', 1.5); 
                txt = 'line node'; 
                text(-n(1)*a*2, -n(2)*a*2, -n(3)*a*2, txt)

            % inclination change 1
            c = [0.4940 0.1840 0.5560]; 
            plot_traj(X_Nsat_ei1, c)
            
            % inclination change 2
            c = [0.4660 0.6740 0.1880]; 
            plot_traj(X_Nsat_ei2, c)
            
            % inclination change 3
            c = [0.3010 0.7450 0.9330]; 
            plot_traj(X_Nsat_ei3, c)

            % inclination change 4
            c = [0.6350 0.0780 0.1840]; 
            plot_traj(X_Nsat_ei4, c)

            title('Raise peri, incl change') 
%         legend('Neptune', 'sat', 'Triton', 'ecliptic', 'J2000_X', 'J2000_Y', 'J2000_Z'); 

        % inclination change, RAAN change 
        subplot(2,2,4) 
            plot_NT(const, X_NT, rnorm_T)
            
            % inclination change 
            c = [0.4940 0.1840 0.5560]; 
            plot_traj(X_Nsat_ei, c)

        sgtitle(ftitle)
    
end

% ------------------------------------------------------------------------
% plot NEPTUNE and TRITON 

function plot_NT(const, X_NT, rnorm_T)

    % Plot Neptune + Triton 
    ellipsoid(0, 0, 0, const.RN + 1000, const.RN + 1000, const.RN + 1000); alpha 0.2; shading interp; 
    hold on; grid on; axis equal 
    ellipsoid(0, 0, 0, const.RN, const.RN, const.RN); alpha 0.4; shading interp; 
    plot3(X_NT(:,1), X_NT(:,2), X_NT(:,3), 'b', 'linewidth', 1.2);
    plot3(X_NT(1,1), X_NT(1,2), X_NT(1,3), 'bo'); 
    plot3(X_NT(end,1), X_NT(end,2), X_NT(end,3), 'b^'); 
%         patch([1 -1 -1 1]*rnorm_T, [1 1 -1 -1]*rnorm_T, [0 0 0 0]*rnorm_T, [1 1 -1 -1]*rnorm_T); alpha 0.2  

    % J2000 axes 
    quiver3(0, 0, 0, rnorm_T, 0, 0)
    quiver3(0, 0, 0, 0, rnorm_T, 0); 
    quiver3(0, 0, 0, 0, 0, rnorm_T); 
        txt = 'J2000_x';
        text(rnorm_T, 0, 0, txt)
        
    view(-54, 64)
    
    xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)'); 

end 

% ------------------------------------------------------------------------
% plot start, middle, and end 

function plot_traj(X_Nsat_e1, c)
    plot3(X_Nsat_e1(1,1), X_Nsat_e1(1,2), X_Nsat_e1(1,3), 'o', 'color', c); 
    plot3(X_Nsat_e1(:,1), X_Nsat_e1(:,2), X_Nsat_e1(:,3), '', 'color', c, 'linewidth', 2); 
    plot3(X_Nsat_e1(end,1), X_Nsat_e1(end,2), X_Nsat_e1(end,3), '^', 'color', c); 
end 

%% subfunctions 

function [t, X_Nsat_ei] = incl_change(X_Nsat_ep, OE_T, dt, const, options)

% obtain OE for input state (orbit)
oe_di = rvOrb.rv2orb(X_Nsat_ep(end,:), const.muN); 

% orbit normal 
h = cross(X_Nsat_ep(end, 1:3), X_Nsat_ep(end, 4:6)); 
h = h / norm(h); 

% desired plane change direction 
h_des = -h; 

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