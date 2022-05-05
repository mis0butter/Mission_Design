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

% % hyperbolic orbit 
% % rv0_test = [19813.3, -16908.2, 2612.7, -22.953, -13.324, 11.316]; 
% 
% % "beginning" of hyperbolic orbit 
% rv0_sat = [ 
%           610470.376325053
%           739754.071895096
%          -468914.276897989
%          -9.72778840140298
%          -12.8560823075242
%           7.93760951868089]; 
% 
% OE0_sat = rvOrb.rv2orb(rv0_sat, const.muN); 
% 
% % satellite period 
% T_sat = 2*pi*sqrt( OE0_sat(1)^3 / const.muN ); 

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
[t, X_Nsat_h] = ode45(@fn.EOM, [0 : -dt : -T0_T/8], rv0_sat, options); 

rv0_sat = X_Nsat_h(end,:); 

%% integrate EOM - satellite 

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
rv_Nsat_min = X_Nsat_h(i_min, :); 

% augment OE "delta v" aerocapture from hyperbolic --> elliptical 
OE_Nsat_min = rvOrb.rv2orb(rv_Nsat_min, const.muN); 
OE_Nsat_min(1) = OE_T(1); 
OE_Nsat_min(2) = 0.9; 

% obtain state of augmented elliptical orbit 
rv0_test = rvOrb.orb2rv(OE_Nsat_min, const.muN); 

% extract velocity direction, save over periapsis state 
rv_Nsat_min(4:6) = rv0_test(4:6)*1.2; 
OE_Nsat_min = rvOrb.rv2orb(rv_Nsat_min, const.muN); 
T_Nsat_min = 2*pi*sqrt( OE_Nsat_min(1)^3 / const.muN ); 

% change RAAN to match Triton's 
% OE_Nsat_min(5) = OE_T(end,5); 
% X_Nsat_min = rvOrb.orb2rv(OE_Nsat_min, const.muN); 

% propagate elliptical orbit 
[t, X_Nsat_e1] = ode45(@fn.EOM, [0: dt : T_Nsat_min], rv_Nsat_min, options); 

% ------------------------------------------------------------------------ 
% 2nd aerocapture ellipse (decrease apoapsis) 

% obtain state of augmented elliptical orbit 
rv0_test = rvOrb.orb2rv(OE_Nsat_min, const.muN); 

% extract velocity direction, save over periapsis state 
rv_Nsat_min(4:6) = rv0_test(4:6)*0.98; 
OE_Nsat_min = rvOrb.rv2orb(rv_Nsat_min, const.muN); 
T_Nsat_min = 2*pi*sqrt( OE_Nsat_min(1)^3 / const.muN ); 

% propagate elliptical orbit 
[t, X_Nsat_e2] = ode45(@fn.EOM, [0: dt : T_Nsat_min*0.9], rv_Nsat_min, options); 

% ------------------------------------------------------------------------ 
% raise periapsis  

% compute rnorm again 
rnorm_Nsat = []; 
for i = 1:length(X_Nsat_e2)
    rnorm_Nsat(i,:) = norm(X_Nsat_e2(i,1:3)); 
end 

% now find index of MAX norm - save state at apoapsis 
i_max = find(rnorm_Nsat == max(rnorm_Nsat)); 
rv_Nsat_max = X_Nsat_e2(i_max, :); 

% augment OE "delta v" raise periapsis 
OE_Nsat_max = rvOrb.rv2orb(rv_Nsat_max, const.muN); 

% it doesn't make sense, but decrease semimajor 
OE_Nsat_max(1) = OE_Nsat_max(1)*0.2;        
% OE_Nsat_max(2) = 0.01; 

% obtain state of raised periapsis elliptical orbit 
rv0_test = rvOrb.orb2rv(OE_Nsat_max, const.muN); 
rv_Nsat_max(4:6) = rv0_test(4:6); 
OE_Nsat_max = rvOrb.rv2orb(rv_Nsat_max, const.muN); 

% keep trying ... 
OE_Nsat_max(2) = 0.01; 
rv0_test = rvOrb.orb2rv(OE_Nsat_max, const.muN); 
rv_Nsat_max(4:6) = rv0_test(4:6); 

OE_Nsat_max = rvOrb.rv2orb(rv_Nsat_max, const.muN); 
OE_Nsat_max(2) = 0.01; 
rv0_test = rvOrb.orb2rv(OE_Nsat_max, const.muN); 
rv_Nsat_max(4:6) = rv0_test(4:6); 

OE_Nsat_max = rvOrb.rv2orb(rv_Nsat_max, const.muN); 
OE_Nsat_max(2) = 0.01; 
rv0_test = rvOrb.orb2rv(OE_Nsat_max, const.muN); 
rv_Nsat_max(4:6) = rv0_test(4:6); 

OE_Nsat_max = rvOrb.rv2orb(rv_Nsat_max, const.muN); 
OE_Nsat_max(2) = 0.01; 
rv0_test = rvOrb.orb2rv(OE_Nsat_max, const.muN); 
rv_Nsat_max(4:6) = rv0_test(4:6); 

OE_Nsat_max = rvOrb.rv2orb(rv_Nsat_max, const.muN); 
T_Nsat_max = 2*pi*sqrt( OE_Nsat_max(1)^3 / const.muN ); 

% propagate raised periapsis elliptical orbit 
[t, X_Nsat_ep] = ode45(@fn.EOM, [0: dt : T_Nsat_max * 0.8], rv_Nsat_max, options); 

%% test lambert 

X_start = rv_Nsat_max'; 

phi_des = 75; 
plot_option = 1; 
[ell_1_min, ell_2_min, amin_AU, emin] = lambert_prob(X_start, t0, phi_des, plot_option);

% % BETTADPUR LAMBERT FN 
% [ell_1, ell_2] = bett_lambert_NT(X_start, t0, phi_des, const, plot_option); 

phi_d_hist  = []; 
phi_a_hist  = []; 
tof_hist = []; 
phi_t_hist = []; 

phi0 = 15; 
plot_option = 0; 
for phi = [ phi0 : 15 : 165] 
    
    phi_t_des = phi; 
    [ell_1_min, ell_2_min, amin_AU, emin] = lambert_prob(X_start, t0, phi_t_des, plot_option);
%     [ell_1_min, ell_2_min] = lambert_prob(t0, phi_t_des, 0); 
    
    % departure angle: 
    % (1) ELL 1 SHORT, (2) ELL 1 LONG, (3) ELL 2 SHORT, (4) ELL 2 LONG 
    phi_d_hist = [phi_d_hist; ... 
        ell_1_min.phi_ds, ell_1_min.phi_dl, ell_2_min.phi_ds, ell_2_min.phi_dl]; 

    % arrival angle 
    % (1) ELL 1 SHORT, (2) ELL 1 LONG, (3) ELL 2 SHORT, (4) ELL 2 LONG 
    phi_a_hist = [phi_a_hist; ... 
        ell_1_min.phi_as, ell_1_min.phi_al, ell_2_min.phi_as, ell_2_min.phi_al]; 

    % time of flight 
    % (1) ELL 1 SHORT, (2) ELL 1 LONG, (3) ELL 2 SHORT, (4) ELL 2 LONG 
    tof_hist = [tof_hist; ... 
        ell_1_min.dt_s, ell_1_min.dt_l, ell_2_min.dt_s, ell_2_min.dt_l]; 

    % transfer angle 
    phi_t_hist = [phi_t_hist; phi_t_des]; 
        
end 


colors = {'k', 'b', 'r', 'g'}; 
style = {'p', '^', '--', '.'}; 
lwidth = [3, 1.5, 1, 1]; 

figure('position', [100 100 700 700])
    subplot(3,1,1) 
        hold on;
        for i = 1:4
%             scatter(phi_t_hist, phi_d_hist(:,i), 4, colors{i});  
            h(i) = plot(phi_t_hist, phi_d_hist(:,i), [colors{i} style{i}]); 
        end 
        title('Departure Angle')
        ylabel('deg') 
        legend(h, 'ell 1 short', 'ell 1 long', 'ell 2 short', 'ell 2 long', 'location', 'eastoutside'); 
    
    subplot(3,1,2) 
        hold on;
        for i = 1:4
%             scatter(phi_t_hist, phi_a_hist(:,i), 4, colors{i});  
            h(i) = plot(phi_t_hist, phi_a_hist(:,i), [colors{i} style{i}]); 
        end 
        title('Arrival Angle') 
        ylabel('deg') 
        legend(h, 'ell 1 short', 'ell 1 long', 'ell 2 short', 'ell 2 long', 'location', 'eastoutside'); 
        
    subplot(3,1,3) 
        hold on; 
        for i = 1:4
%             scatter(phi_t_hist, tof_hist(:,i), 4, colors{i});  
            h(i) = plot(phi_t_hist, tof_hist(:,i), [colors{i} style{i}]); 
        end 
        title('Time of Flight') 
        ylabel('days') 
        legend(h, 'ell 1 short', 'ell 1 long', 'ell 2 short', 'ell 2 long', 'location', 'eastoutside'); 
        
    xlabel('Transfer Angle Phi (deg)') 
    sgtitle('Min Energy Parameters') 

%%

% ------------------------------------------------------------------------ 
% 1st RAAN change maneuver 

% obtain desired (Triton) RAAN 
oe_ep = rvOrb.rv2orb(X_Nsat_ep(end,:), const.muN); 
oe_des = oe_ep; 
oe_des(5) = OE_T(end,5); 

% convert desired OEs back to state 
X_Nsat_er = rvOrb.orb2rv(oe_des, const.muN); 

% propagate desired changed RAAN elliptical orbit 
[t, X_Nsat_er] = ode45(@fn.EOM, [0: dt : T_Nsat_max * 0.8], X_Nsat_er, options); 

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

% [t, X_Nsat_RAAN] = RAAN_change(X_Nsat_ep, const, OE_T, dt, options); 

% ------------------------------------------------------------------------ 
% 1st inclination change maneuver 

[t, X_Nsat_ei1] = incl_change(X_Nsat_RAAN, OE_T, dt, const, options); 
[t, X_Nsat_ei2] = incl_change(X_Nsat_ei1, OE_T, dt, const, options); 
[t, X_Nsat_ei3] = incl_change(X_Nsat_ei2, OE_T, dt, const, options); 
[t, X_Nsat_ei4] = incl_change(X_Nsat_ei3, OE_T, dt, const, options); 

%% test lambert again 

X_start = X_Nsat_ei4(round(end/2),:)'; 

phi_des = 75; 
plot_option = 1; 
[ell_1_min, ell_2_min, amin_AU, emin] = lambert_prob(X_start, t0, phi_des, plot_option);
%%
% % BETTADPUR LAMBERT FN 
% [ell_1, ell_2] = bett_lambert_NT(X_start, t0, phi_des, const, plot_option); 

phi_d_hist  = []; 
phi_a_hist  = []; 
tof_hist = []; 
phi_t_hist = []; 

phi0 = 15; 
plot_option = 0; 
for phi = [ phi0 : 15 : 165] 
    
    phi_t_des = phi; 
    [ell_1_min, ell_2_min, amin_AU, emin] = lambert_prob(X_start, t0, phi_t_des, plot_option);
%     [ell_1_min, ell_2_min] = lambert_prob(t0, phi_t_des, 0); 
    
    % departure angle: 
    % (1) ELL 1 SHORT, (2) ELL 1 LONG, (3) ELL 2 SHORT, (4) ELL 2 LONG 
    phi_d_hist = [phi_d_hist; ... 
        ell_1_min.phi_ds, ell_1_min.phi_dl, ell_2_min.phi_ds, ell_2_min.phi_dl]; 

    % arrival angle 
    % (1) ELL 1 SHORT, (2) ELL 1 LONG, (3) ELL 2 SHORT, (4) ELL 2 LONG 
    phi_a_hist = [phi_a_hist; ... 
        ell_1_min.phi_as, ell_1_min.phi_al, ell_2_min.phi_as, ell_2_min.phi_al]; 

    % time of flight 
    % (1) ELL 1 SHORT, (2) ELL 1 LONG, (3) ELL 2 SHORT, (4) ELL 2 LONG 
    tof_hist = [tof_hist; ... 
        ell_1_min.dt_s, ell_1_min.dt_l, ell_2_min.dt_s, ell_2_min.dt_l]; 

    % transfer angle 
    phi_t_hist = [phi_t_hist; phi_t_des]; 
        
end 


colors = {'k', 'b', 'r', 'g'}; 
style = {'p', '^', '--', '.'}; 
lwidth = [3, 1.5, 1, 1]; 

figure('position', [100 100 700 700])
    subplot(3,1,1) 
        hold on;
        for i = 1:4
%             scatter(phi_t_hist, phi_d_hist(:,i), 4, colors{i});  
            h(i) = plot(phi_t_hist, phi_d_hist(:,i), [colors{i} style{i}]); 
        end 
        title('Departure Angle')
        ylabel('deg') 
        legend(h, 'ell 1 short', 'ell 1 long', 'ell 2 short', 'ell 2 long', 'location', 'eastoutside'); 
    
    subplot(3,1,2) 
        hold on;
        for i = 1:4
%             scatter(phi_t_hist, phi_a_hist(:,i), 4, colors{i});  
            h(i) = plot(phi_t_hist, phi_a_hist(:,i), [colors{i} style{i}]); 
        end 
        title('Arrival Angle') 
        ylabel('deg') 
        legend(h, 'ell 1 short', 'ell 1 long', 'ell 2 short', 'ell 2 long', 'location', 'eastoutside'); 
        
    subplot(3,1,3) 
        hold on; 
        for i = 1:4
%             scatter(phi_t_hist, tof_hist(:,i), 4, colors{i});  
            h(i) = plot(phi_t_hist, tof_hist(:,i), [colors{i} style{i}]); 
        end 
        title('Time of Flight') 
        ylabel('days') 
        legend(h, 'ell 1 short', 'ell 1 long', 'ell 2 short', 'ell 2 long', 'location', 'eastoutside'); 
        
    xlabel('Transfer Angle Phi (deg)') 
    sgtitle('Min Energy Parameters') 

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

    ftitle = 'Hyperb, 1st aero, 2nd aero'; 
    figure('name', ftitle, 'position', pos); 
    
        plot_NTsun(const, X_NT, rnorm_T, r_Nsun)

        % hyperbolic trajectory 
        c = [1 0 1]; 
        plot_traj(X_Nsat_h, c, 'hyperb'); 
        plot3(rv_Nsat_min(1), rv_Nsat_min(2), rv_Nsat_min(3), 'mp'); 

        % 1st elliptical aerocapture 
        c = [0 1 1]; 
        plot_traj(X_Nsat_e1, c, 'aero'); 

        % 2nd elliptical aerocapture 
        plot_traj(X_Nsat_e2, c); 

        % APOAPSIS 
        plot3(X_Nsat_e2(i_max,1), X_Nsat_e2(i_max,2), X_Nsat_e2(i_max,3), 'p'); 

        title(ftitle)

% ------------------------------------------------------------------------
    % 2nd aerocapture, raise periapsis 

    ftitle = '2nd aero, raise peri'; 
    pos = pos + [25 25 0 0]; 
    figure('name', ftitle, 'position', pos); 
        plot_NTsun(const, X_NT, rnorm_T, r_Nsun)

        % 2nd elliptical aerocapture 
        c = [0 1 1]; 
        plot_traj(X_Nsat_e2, c, 'aero'); 

        % APOAPSIS 
        plot3(X_Nsat_e2(i_max,1), X_Nsat_e2(i_max,2), X_Nsat_e2(i_max,3), 'p'); 

        % raise periapsis 
        c = [0.4660 0.6740 0.1880]; 
        plot_traj(X_Nsat_ep, c, 'raise peri')

        title(ftitle) 

% ------------------------------------------------------------------------
    % raise periapsis, RAAN change 

    ftitle = 'Raise peri, RAAN change'; 
    pos = pos + [25 25 0 0]; 
    figure('name', ftitle, 'position', pos); 

        plot_NTsun(const, X_NT, rnorm_T, r_Nsun)
        c = [0.4660 0.6740 0.1880]; 
        plot_traj(X_Nsat_ep, c, 'raise peri')
        % plot_traj(X_Nsat_cur, [0 1 1], 'curr circ'); 

        a = OE_T(1); 
        n = line_node(X_Nsat_ep) * a * 2; 
        quiver3(0, 0, 0, n(1), n(2), n(3), 'k', 'linewidth', 1.5); 
            txt = 'Curr asc. node'; 
            text(n(1), n(2), n(3), txt)

        % plot_traj(X_Nsat_des, [1 0 1], 'des circ'); 
        % plot_traj(X_Nsat_er, [1 0 1], 'des ell')

        n = line_node(X_NT) * a * 2; 
        quiver3(0, 0, 0, n(1), n(2), n(3), 'k', 'linewidth', 1.5); 
            txt = 'Des asc. node'; 
            text(n(1), n(2), n(3), txt)

        % state at burn (initial) 
        quiver3(rv_raan_i(1), rv_raan_i(2), rv_raan_i(3), ... 
            v_raan_i(1)*30000, ... 
            v_raan_i(2)*30000, ... 
            v_raan_i(3)*30000, 'r'); 
            txt = 'vel init'; 
            text(rv_raan_i(1) + v_raan_i(1)*30000, ... 
                rv_raan_i(2) + v_raan_i(2)*30000, ... 
                rv_raan_i(3) + v_raan_i(3)*30000, txt)  

        % final orbit 
        rv_raan_f = rvOrb.orb2rv(oe_raan_f, const.muN); 
        r_raan_f = rv_raan_f(1:3) / norm(rv_raan_f(1:3)) * OE_T(1) * 1.4; 
        v_raan_f = rv_raan_f(4:6); 
        quiver3(0, 0, 0, r_raan_f(1), r_raan_f(2), r_raan_f(3), 'r', 'linewidth', 1.5)
            txt = 'Burn point'; 
            text(r_raan_f(1), r_raan_f(2), r_raan_f(3), txt)
        quiver3(rv_raan_i(1), rv_raan_i(2), rv_raan_i(3), ... 
            v_raan_f(1)*30000, ... 
            v_raan_f(2)*30000, ... 
            v_raan_f(3)*30000, 'r'); 
            txt = 'vel fin'; 
            text(rv_raan_f(1) + v_raan_f(1)*30000, ... 
                rv_raan_f(2) + v_raan_f(2)*30000, ... 
                rv_raan_f(3) + v_raan_f(3)*30000, txt)   

        % RAAN maneuver 
        c = [0.6350 0.0780 0.1840]; 
        plot_traj(X_Nsat_RAAN, c, 'change RAAN')

        title(ftitle) 

% ------------------------------------------------------------------------
    % RAAN change, inclination change

    ftitle = 'RAAN change, incl change'; 
    pos = pos + [25 25 0 0]; 
    figure('name', ftitle, 'position', pos); 
    
        plot_NTsun(const, X_NT, rnorm_T, r_Nsun)
        
        % RAAN change 
        c = [0.6350 0.0780 0.1840]; 
        plot_traj(X_Nsat_RAAN, c, 'change RAAN'); 

        % inclination change 1
        c = [0.3010 0.7450 0.9330]; 
        plot_traj(X_Nsat_ei1, c, 'change incl')

        % inclination change 2
        plot_traj(X_Nsat_ei2, c)

        % inclination change 3
        plot_traj(X_Nsat_ei3, c)

        % inclination change 4
        plot_traj(X_Nsat_ei4, c)
        
        n = line_node(X_Nsat_ei1) * a * 2; 
        quiver3(0, 0, 0, n(1), n(2), n(3), 'k', 'linewidth', 1.5); 
            txt = 'asc. node'; 
            text(n(1), n(2), n(3), txt)

        title(ftitle)
    
end

% ------------------------------------------------------------------------
% plot NEPTUNE and TRITON 

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
    plot3(X_NT(:,1), X_NT(:,2), X_NT(:,3), 'b', 'linewidth', 1);
    lgd_T = plot3(X_NT(1:end/6,1), X_NT(1:end/6,2), X_NT(1:end/6,3), 'b', 'linewidth', 2);
    plot3(X_NT(1,1), X_NT(1,2), X_NT(1,3), 'bo'); 
    plot3(X_NT(end/6,1), X_NT(end/6,2), X_NT(end/6,3), 'b^'); 

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

function plot_traj(X_Nsat_e1, c, leg_str)
    
    if nargin < 3
        plot3(X_Nsat_e1(:,1), X_Nsat_e1(:,2), X_Nsat_e1(:,3), '', 'color', c, 'linewidth', 2); 
    else
        hleg = get(gca, 'Legend'); 
        hleg.AutoUpdate = 'on'; 
        plot3(X_Nsat_e1(:,1), X_Nsat_e1(:,2), X_Nsat_e1(:,3), '', 'color', c, 'linewidth', 2); 
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