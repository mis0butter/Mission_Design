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
for i = 0 : dt : T0_T
    
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

target = 'Sun'; 
X_Nsun  = spice_state(et, target, frame, abcorr, observer); 

% sun vector with Triton semimajor axis length 
r_Nsun = X_Nsun(1:3) / norm(X_Nsun(1:3)) * OE_T(1); 

%% initial state 

% desired OEs at perigee 
OE_Nsat_min = [ 
    153727.013558483
    0.83852878211854
    2.67992509051225
    0.618907783952636
    1.5
    6.27786334026983 ]; 
      
rv0_sat = rvOrb.orb2rv(OE_Nsat_min, const.muN); 
  
% change to hyperbolic 
OE_h = OE_Nsat_min; 
OE_h(1) = -22045.6960334535; 
OE_h(2) = 2.12595341736863; 

rv0_sat = rvOrb.orb2rv(OE_h, const.muN)

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% propagate negative hyperbolic 
[t, X_Nsat_h] = ode45(@fn.EOM, [0 : -dt : -T0_T/20], rv0_sat, options); 

rv0_sat = X_Nsat_h(end,:); 
X_Nsat_h = flip(X_Nsat_h); 

%% integrate EOM - satellite 

% Set run state 
disp('Running sim ...')

dt = 10; 

% propagate hyperbolic orbit 
tic
[t, X_Nsat_h_orbit] = ode45(@fn.EOM, [0 : dt : T0_T/10], rv0_sat, options); 
toc 
disp('Pos and Vel end: ')
disp(X_Nsat_h_orbit(end, 1:6)'); 

% ------------------------------------------------------------------------ 
% 1st aerocapture ellipse (decrease apoapsis) 

% compute rnorm 
rnorm_Nsat = []; 
for i = 1:length(X_Nsat_h_orbit)
    rnorm_Nsat(i,:) = norm(X_Nsat_h_orbit(i,1:3)); 
end 

% find index of min norm - save state at periapsis 
i_min = find(rnorm_Nsat == min(rnorm_Nsat)); 
rvf_Nsat_h = X_Nsat_h_orbit(i_min, :); 
rvf_Nsat_h = X_Nsat_h(end,:); 
rv0_Nsat_aero1 = rvf_Nsat_h; 

% augment OE "delta v" aerocapture from hyperbolic --> elliptical 
OE_Nsat_min = rvOrb.rv2orb(rv0_Nsat_aero1, const.muN); 
OE_Nsat_min(1) = OE_T(1); 
OE_Nsat_min(2) = 0.9; 

% obtain state of augmented elliptical orbit 
rv0_test = rvOrb.orb2rv(OE_Nsat_min, const.muN); 

% extract velocity direction, save over periapsis state 
rv0_Nsat_aero1(4:6) = rv0_test(4:6)*1.2; 
OE_Nsat_min = rvOrb.rv2orb(rv0_Nsat_aero1, const.muN); 
T_Nsat_min = 2*pi*sqrt( OE_Nsat_min(1)^3 / const.muN ); 

% change RAAN to match Triton's 
% OE_Nsat_min(5) = OE_T(end,5); 
% X_Nsat_min = rvOrb.orb2rv(OE_Nsat_min, const.muN); 

% delta v 
dv_h_aero1 = abs(norm(rvf_Nsat_h(4:6)) - norm(rv0_Nsat_aero1(4:6))); 

% propagate elliptical orbit 
[t, X_Nsat_aero1_orbit] = ode45(@fn.EOM, [0: dt : T_Nsat_min], rv0_Nsat_aero1, options); 


%% test lambert again 

% find apoapsis 
for i = 1:length(X_Nsat_aero1_orbit)
    rnorm(i,:) = norm(X_Nsat_aero1_orbit(i,1:3)); 
end
i_max = find(rnorm == max(rnorm));
rvf_Nsat_aero1 = X_Nsat_aero1_orbit(i_max,:)'; 
X_Nsat_aero1 = X_Nsat_aero1_orbit(1:i_max, :); 
rv0_Nsat_lambert = rvf_Nsat_aero1; 

phi_des = 120; 
plot_option = 1; 
[ell_1_min, ell_2_min, amin_AU, emin] = lambert_prob(rv0_Nsat_lambert, t0, phi_des, plot_option);
X_Nsat_lambert = ell_1_min.rv_s; 
rv0_Nsat_lambert = X_Nsat_lambert(1,:); 

% delta v 
dv_aero1_lambert = abs(norm(rvf_Nsat_aero1(4:6)) - norm(rv0_Nsat_lambert(4:6))); 


%% transition into science orbit 

rvf_Nsat_lambert = X_Nsat_lambert(end,:); 
rv0_Nsat_sci = rvf_Nsat_lambert; 

% find same point in Triton orbit 
% calculate angle between Triton vector and final lambert vector 
for i = 1:length(X_NT)
    
    a = X_NT(i,1:3); 
    b = rv0_Nsat_sci(1:3); 
    phi(i,:) = acosd( dot( a, b)/(norm(a)*norm(b)) ); 
    
end 
i_min = find(phi == min(phi)); 

% go into science orbit 
rv0_Nsat_sci(4:6) = X_NT(i_min, 4:6); 
rv0_Nsat_sci(4:6) = rv0_Nsat_sci(4:6)*0.9; 
OE_Nsat_sci = rvOrb.rv2orb(rv0_Nsat_sci, const.muN); 
T_Nsat_sci = 2*pi*sqrt( OE_Nsat_sci(1)^3 / const.muN ); 

% delta v 
dv_lambert_sci = norm(rvf_Nsat_lambert(4:6)) - norm(rv0_Nsat_sci(4:6)); 

% propagate elliptical orbit 
[t, X_Nsat_sci] = ode45(@fn.EOM, [0: dt : T_Nsat_sci], rv0_Nsat_sci, options); 

clear rnorm 
for i = 1:length(X_Nsat_sci)
    rnorm(i,:) = norm(X_Nsat_sci(i, 1:3)); 
end 
i_min = find(rnorm == min(rnorm)); 

while abs(min(rnorm) - const.RN) > 1000

    % go into science orbit 
    if abs(min(rnorm) - const.RN) > 10000 
        SF = 0.98; 
    else
        SF = 0.99; 
    end 
    rv0_Nsat_sci(4:6) = rv0_Nsat_sci(4:6)*SF;
    OE_Nsat_sci = rvOrb.rv2orb(rv0_Nsat_sci, const.muN); 
    T_Nsat_sci = 2*pi*sqrt( OE_Nsat_sci(1)^3 / const.muN ); 

    % delta v 
    dv_lambert_sci = norm(rvf_Nsat_lambert(4:6)) - norm(rv0_Nsat_sci(4:6)); 

    % propagate elliptical orbit 
    [t, X_Nsat_sci] = ode45(@fn.EOM, [0: dt : T_Nsat_sci], rv0_Nsat_sci, options); 

    clear rnorm 
    for i = 1:length(X_Nsat_sci)
        rnorm(i,:) = norm(X_Nsat_sci(i, 1:3)); 
    end 
    i_min = find(rnorm == min(rnorm)); 
    abs(min(rnorm) - const.RN)
    
end 

%% plot 

% Neptune sphere 
[X, Y, Z] = sphere; 
XN = X * const.RN; YN = Y * const.RN; ZN = Z * const.RN; 

% axis length 
rnorm_T = norm(X_NT(end,1:3)); 

pos = [100 100 800 600]; 
plot_orbit = 1; 
if plot_orbit == 1

% ------------------------------------------------------------------------
    % hyperbolic, 1st aerocapture, 2nd aerocapture 
    
    ftitle1 = 'Hyperb --> Aerocapture'; 
    ftitle2 = sprintf('dv = %.6g km/s', dv_h_aero1); 
    figure('name', ftitle1, 'position', pos); 
        % plot Neptune, sun 
        plot_NTsun(const, X_NT, rnorm_T, r_Nsun)
        % hyperbolic trajectory 
        c = [1 0 1]; 
        plot_traj(X_Nsat_h, c, 'hyperb', '-'); 
        plot_traj(X_Nsat_h_orbit, c); 
        % 1st elliptical aerocapture 
        c = [0 1 1]; 
        plot_traj(X_Nsat_aero1_orbit, c); 
        plot_traj(X_Nsat_aero1, c, 'aero', '-'); 
        % title 
        title( {ftitle1; ftitle2} )
        
    ftitle1 = 'Aerocapture --> Lambert'; 
    ftitle2 = sprintf('dv = %.6g km/s', dv_aero1_lambert); 
    figure('name', ftitle1, 'position', pos + [0 0 100 0]); 
        % plot Neptune, sun 
        plot_NTsun(const, X_NT, rnorm_T, r_Nsun)
        % 1st elliptical aerocapture 
        c = [0 1 1]; 
        plot_traj(X_Nsat_aero1_orbit, c); 
        plot_traj(X_Nsat_aero1, c, 'aero', '-'); 
        % plot lambert 
        c = [1 0 0]; 
        plot_traj(X_Nsat_lambert, c, 'lambert', '-'); 
        % title 
        title( {ftitle1; ftitle2} )
            
    ftitle1 = 'Lambert --> Science'; 
    ftitle2 = sprintf('dv = %.6g km/s', dv_lambert_sci); 
    figure('name', ftitle1, 'position', pos + [0 0 100 0]);
        % plot Neptune, sun 
        plot_NTsun(const, X_NT, rnorm_T, r_Nsun)
        % plot lambert 
        c = [1 0 0]; 
        plot_traj(X_Nsat_lambert, c, 'lambert', '-'); 
        % plot science orbit
        c = [0.4660 0.6740 0.1880]; 
        plot_traj(X_Nsat_sci, c, 'science', '-'); 
        % title 
        title( {ftitle1; ftitle2} )
    
end


%% plot NEPTUNE and TRITON 

function plot_NTsun(const, X_NT, rnorm_T, X_Nsun)

    % Neptune sphere 
    [X, Y, Z] = sphere; 
    XN = X * const.RN; YN = Y * const.RN; ZN = Z * const.RN; 
    XN_atm = X * (const.RN + 1000); YN_atm = Y * (const.RN + 1000); ZN_atm = Z * (const.RN + 1000); 

    % Plot Neptune 
    lgd_n = surf(XN, YN, ZN); alpha 0.2; shading interp; 
    hold on; grid on; axis equal 
    surf(XN_atm, YN_atm, ZN_atm); alpha 0.2; shading interp; 
    
    % Plot Triton 
    lgd_T = plot3(X_NT(:,1), X_NT(:,2), X_NT(:,3), 'b--', 'linewidth', 2);
%     lgd_T = plot3(X_NT(1:end/6,1), X_NT(1:end/6,2), X_NT(1:end/6,3), 'b--', 'linewidth', 2);
%     plot3(X_NT(1,1), X_NT(1,2), X_NT(1,3), 'bo'); 
%     plot3(X_NT(end/6,1), X_NT(end/6,2), X_NT(end/6,3), 'b^'); 

    % axes SF 
    SF = 1.4; 
    
    % J2000 axes 
    rnorm_T = rnorm_T*SF; 
    quiver3(0, 0, 0, rnorm_T, 0, 0)
    quiver3(0, 0, 0, 0, rnorm_T, 0); 
    quiver3(0, 0, 0, 0, 0, rnorm_T); 
        txt = 'x_{ecl}';
        text(rnorm_T, 0, 0, txt)
        
    % SUN 
    c = [0.9290 0.6940 0.1250]; 
    X_Nsun = X_Nsun*SF; 
    lgd_S = quiver3(0, 0, 0, X_Nsun(1), X_Nsun(2), X_Nsun(3), 'color', c); 
        
    view(-54, 64)
    
%     legend('Neptune', 
    h_lgd = [lgd_n(1), lgd_T(1), lgd_S(1)]; 
    legend(h_lgd, 'Neptune', 'Triton', 'Sun', 'AutoUpdate', 'off'); 
    xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)'); 

end 

% ------------------------------------------------------------------------
% plot start, middle, and end 

function plot_traj(X_Nsat_e1, c, leg_str, style)
    
    if ~exist('style', 'var')
        style = '--'; 
    end 
    if ~exist('leg_str', 'var')
        plot3(X_Nsat_e1(:,1), X_Nsat_e1(:,2), X_Nsat_e1(:,3), style, 'color', c, 'linewidth', 2); 
    else
        hleg = get(gca, 'Legend'); 
        hleg.AutoUpdate = 'on'; 
        plot3(X_Nsat_e1(:,1), X_Nsat_e1(:,2), X_Nsat_e1(:,3), style, 'color', c, 'linewidth', 2); 
        i = length(hleg.String); 
        hleg.String{i} = leg_str; 
        hleg.AutoUpdate = 'off'; 
    end 
    
    
    plot3(X_Nsat_e1(1,1), X_Nsat_e1(1,2), X_Nsat_e1(1,3), 'o', 'color', c); 
    
    plot3(X_Nsat_e1(end,1), X_Nsat_e1(end,2), X_Nsat_e1(end,3), '^', 'color', c); 
    
    
end 

% obtain line node 
function n = line_node(X_Nsat_ep)

    % orbit normal 
    h = cross(X_Nsat_ep(end, 1:3), X_Nsat_ep(end, 4:6)); 
    h = h / norm(h); 

    % line node 
    n = cross([0 0 1], h); 
    n = n / norm(n); 

end 

%% subfunctions 

function [t, X_Nsat_RAAN] = RAAN_change(X_Nsat_ep, const, OE_T, dt, options)

% obtain desired (Triton) RAAN 
oe_ep = rvOrb.rv2orb(X_Nsat_ep(end,:), const.muN); 
oe_des = oe_ep; 
oe_des(5) = OE_T(end,5); 

% convert desired OEs back to state 
X_Nsat_er = rvOrb.orb2rv(oe_des, const.muN); 

% propagate desired changed RAAN elliptical orbit 
[t, X_Nsat_er] = ode45(@fn.EOM, [0: dt : oe_ep(1) * 0.8], X_Nsat_er, options); 

% convert to (more) circular 
oe_cur = oe_ep;
oe_cur(2) = 0.0001; 
% X_Nsat_cur = rvOrb.orb2rv(oe_cur, const.muN); 
% [t, X_Nsat_cur] = ode45(@fn.EOM, [0: dt : T_Nsat_max*0.9], X_Nsat_cur, options); 

oe_des(2) = 0.0001; 
% X_Nsat_des = rvOrb.orb2rv(oe_des, const.muN); 
% [t, X_Nsat_des] = ode45(@fn.EOM, [0: dt : T_Nsat_max*0.9], X_Nsat_des, options); 

% delta RAAN, incl, arg peri 
d_RAAN = oe_cur(5) - oe_des(5); 
i = oe_cur(3); 

% burn angle 
vnu = acos( cos(i)*cos(i) + sin(i)*sin(i) * cos(d_RAAN) ); 

% solve for location of burn (initial) 
u_i = acos( tan(i) * ( cos(d_RAAN) - cos(vnu) )/sin(vnu) ); 
w = oe_cur(4); 
nu_i = u_i - w; 

% burn at nu (initial) 
oe_raan_i = oe_ep; 
oe_raan_i(6) = nu_i ; 
rv_raan_i = rvOrb.orb2rv(oe_raan_i, const.muN); 
r_raan_i = rv_raan_i(1:3); 
v_raan_i = rv_raan_i(4:6); 

% solve for location of burn (final)
u_f = acos( cos(i) * sin(i) * ( 1 - cos(d_RAAN) )/sin(vnu) ); 
w = oe_des(4); 
nu_f = u_f - w; 

% burn at nu (final) 
oe_raan_f = oe_des; 
oe_raan_f(6) = nu_f ;    
rv_raan_f = rvOrb.orb2rv(oe_raan_f, const.muN); 
r_raan_f = rv_raan_f(1:3); 
v_raan_f = rv_raan_f(4:6); 

% delta v 
dv_raan = v_raan_f - v_raan_i; 
dv_raan_mag = 2*norm(v_raan_i) * sin(vnu/2); 
dv_raan = dv_raan / norm(dv_raan) * dv_raan_mag; 
rv_raan_i(4:6) = rv_raan_i(4:6) + dv_raan; 

oe_raan_i = rvOrb.rv2orb(rv_raan_i, const.muN); 
T_Nsat_raan = 2*pi*sqrt( oe_raan_i(1)^3 / const.muN ); 

% propagate RAAN burn 
[t, X_Nsat_RAAN] = ode45(@fn.EOM, [0: dt : T_Nsat_raan*0.9], rv_raan_i, options); 

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