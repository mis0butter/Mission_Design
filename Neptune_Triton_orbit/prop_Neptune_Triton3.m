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
T_Nsat_max = 2*pi*sqrt( OE_Nsat_max(1)^3 / const.muN ); 

% obtain state of raised periapsis elliptical orbit 
rv0_test = rvOrb.orb2rv(OE_Nsat_max, const.muN); 
X_Nsat_max(4:6) = rv0_test(4:6); 

% propagate raised periapsis elliptical orbit 
[t, X_Nsat_ep] = ode45(@fn.EOM, [0: dt : T_Nsat_max * 3], X_Nsat_max, options); 

% ------------------------------------------------------------------------ 
% inclination change maneuver 

% obtain OEs 
oe_di = rvOrb.rv2orb(X_Nsat_ep(end,:), const.muN); 
w = oe_di(4); 
oe_di(4) = 2*pi - w; 

% orbit normal 
h = cross(X_Nsat_ep(end, 1:3), X_Nsat_ep(end, 4:6)); 
h = h / norm(h); 

% desired plane change direction 
h_des = -h; 

% desired change in inclination = 180 deg - Triton inclination 
di = pi - OE_T(end,3); 

% % line node 
% n = cross([0 0 1], h); 
% n = n / norm(n); 
% 
% k = 0; 
% phi = 1; 
% % find where orbit and line node crossing direction overlap 
% while phi > 0.1 
%     
%     k = k + 1; 
%     r_unit = X_Nsat_ep(k, 1:3) / norm(X_Nsat_ep(k, 1:3));
%     phi    = acos( dot(n, r_unit) ); 
%     
% end 
rv_di = rvOrb.orb2rv(oe_di, const.muN); 

% satellite state and OEs at line of node 
% rv_di = X_Nsat_ep(k,:); 
% oe_di = rvOrb.rv2orb(rv_di, const.muN); 

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
dvi_norm = 1; 

% desired delta velocity vector 
dvi = h_des * dvi_norm; 

% new state 
oe_di(6) = nu; 
rv_di = rvOrb.orb2rv(oe_di, const.muN); 
rv_di(4:6) = rv_di(4:6) + dvi'; 

% fpa again 
fpa = atan( e*sin(nu) / ( 1 + e*cos(nu) ) ); 

% propagate 
[t, X_Nsat_ei] = ode45(@fn.EOM, [0: dt : T_Nsat_min*10], rv_di, options); 

oe_di2 = rvOrb.rv2orb(X_Nsat_ei(end,:), const.muN); 

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
    
        % Plot Neptune + Triton 
        ellipsoid(0, 0, 0, const.RN + 1000, const.RN + 1000, const.RN + 1000); alpha 0.2; shading interp; 
        hold on; grid on; axis equal 
        ellipsoid(0, 0, 0, const.RN, const.RN, const.RN); alpha 0.4; shading interp; 
        plot3(X_NT(:,1), X_NT(:,2), X_NT(:,3), 'b', 'linewidth', 2);
        plot3(X_NT(1,1), X_NT(1,2), X_NT(1,3), 'bo'); 
        plot3(X_NT(end,1), X_NT(end,2), X_NT(end,3), 'b^'); 
%         patch([1 -1 -1 1]*rnorm_T, [1 1 -1 -1]*rnorm_T, [0 0 0 0]*rnorm_T, [1 1 -1 -1]*rnorm_T); alpha 0.2  

        % J2000 axes 
        quiver3(0, 0, 0, v_eq(1), v_eq(2), v_eq(3))
        quiver3(0, 0, 0, 0, rnorm_T, 0); 
        quiver3(0, 0, 0, 0, 0, rnorm_T); 
            txt = 'J2000_x';
            text(v_eq(1),v_eq(2),v_eq(3), txt)
        
        % hyperbolic trajectory 
        plot3(X_Nsat_h(1,1), X_Nsat_h(1,2), X_Nsat_h(1,3), 'mo'); 
        plot3(X_Nsat_h(:,1), X_Nsat_h(:,2), X_Nsat_h(:,3), 'm', 'linewidth', 2); 
        plot3(X_Nsat_h(end,1), X_Nsat_h(end,2), X_Nsat_h(end,3), 'm^'); 
        plot3(X_Nsat_min(1), X_Nsat_min(2), X_Nsat_min(3), 'mp'); 

%         plot3(19813.3, -16908.2, 2612.7, 'p')
        % 1st elliptical aerocapture 
        plot3(X_Nsat_e1(1,1), X_Nsat_e1(1,2), X_Nsat_e1(1,3), 'co'); 
        plot3(X_Nsat_e1(:,1), X_Nsat_e1(:,2), X_Nsat_e1(:,3), 'c', 'linewidth', 2); 
        plot3(X_Nsat_e1(end,1), X_Nsat_e1(end,2), X_Nsat_e1(end,3), 'c^'); 
        
        % 2nd elliptical aerocapture 
        plot3(X_Nsat_e2(1,1), X_Nsat_e2(1,2), X_Nsat_e2(1,3), 'co'); 
        plot3(X_Nsat_e2(:,1), X_Nsat_e2(:,2), X_Nsat_e2(:,3), 'c', 'linewidth', 2); 
        plot3(X_Nsat_e2(end,1), X_Nsat_e2(end,2), X_Nsat_e2(end,3), 'c^'); 

        % APOAPSIS 
        plot3(X_Nsat_e2(i_max,1), X_Nsat_e2(i_max,2), X_Nsat_e2(i_max,3), 'p'); 
        
        % RAAN line node 
        a = OE_Nsat_min(1); 
        quiver3(0, 0, 0, n(1)*a, n(2)*a, n(3)*a, 'k', 'linewidth', 1.5); 
        plot3(rv_di(1), rv_di(2), rv_di(3), 'kp')
            txt = 'RAAN line node'; 
            text(n(1)*a, n(2)*a, n(3)*a, txt)
        quiver3(0, 0, 0, -n(1)*a*2, -n(2)*a*2, -n(3)*a*2, 'k', 'linewidth', 1.5); 
            txt = 'Desc. line node'; 
            text(-n(1)*a*2, -n(2)*a*2, -n(3)*a*2, txt)
        
        % raise periapsis 
        c = [0.8500 0.3250 0.0980]; 
        plot3(X_Nsat_ep(1,1), X_Nsat_ep(1,2), X_Nsat_ep(1,3), 'o', 'color', c); 
        plot3(X_Nsat_ep(:,1), X_Nsat_ep(:,2), X_Nsat_ep(:,3), 'color', c, 'linewidth', 2); 
        plot3(X_Nsat_ep(end,1), X_Nsat_ep(end,2), X_Nsat_ep(end,3), '^', 'color', c); 
        
        % inclination change 
%         plot3(X_Nsat_ep(1,1), X_Nsat_ep(1,2), X_Nsat_ep(1,3), 'go'); 
%         plot3(X_Nsat_ep(:,1), X_Nsat_ep(:,2), X_Nsat_ep(:,3), 'g--', 'linewidth', 2); 
%         plot3(X_Nsat_ep(end,1), X_Nsat_ep(end,2), X_Nsat_ep(end,3), 'g^'); 
        


%         legend('Neptune', 'sat', 'Triton', 'ecliptic', 'J2000_X', 'J2000_Y', 'J2000_Z'); 

        xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)'); 
        title(ftitle)
        
%         view(0, 90); 
        
%     labels = {'a', 'e', 'i', 'w (arg of perigee)', 'O (RAAN)', 'nu (true anomaly)'}; 
%     ftitle = 'Orbital Elements'; 
%     figure('name', ftitle) 
%         for i = 1:6 
%             subplot(6,1,i) 
%             plot(t, OE(:,i))
%             title(labels{i})
%         end 
%         xlabel('time (s)') 
    
end

