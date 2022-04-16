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

%% Obtain initial and final positions 

%  Define parameters for a state lookup:
% t0      = 'Oct 20, 2020 11:00 AM CST'; 
t0      = 'May 22, 2020'; 

%  Convert the epoch to ephemeris time (secs) 
et_t0   = cspice_str2et( t0 );

abcorr  = 'NONE';

% get states --> Sun to Earth 
target   = 'Earth';
frame    = 'ECLIPJ2000'; 
observer = 'Sun';
abcorr   = 'NONE'; 

% get Sun-earth position 
X_sunE  = spice_state(et_t0, target, frame, abcorr, observer); 

% get sun-Mars position 
target   = 'Mars';
X_sunM  = spice_state(et_t0, target, frame, abcorr, observer); 

% get angle between Earth and Mars velocities 
r_sunE = X_sunE(1:3); 
r_sunM = X_sunM(1:3); 
v_sunE = X_sunE(4:6); 
v_sunM = X_sunM(4:6); 

% position and velocity angles 
phi_r = acosd( dot(r_sunE, r_sunM) / (norm(r_sunE)*norm(r_sunM)) ); 
phi_v = acosd( dot(v_sunE, v_sunM) / (norm(v_sunE)*norm(v_sunM)) ); 

% desired phase angle 
phi_des = 179.6234; 

i = 0; 
X_sunM_hist = [X_sunM]; 
X_sunE_hist = [X_sunE]; 
et_hist = [et_t0]; 
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
    
    % get sun-Mars state 
    target = 'Mars';
    X_sunM = spice_state(et, target, frame, abcorr, observer); 
    r_sunM = X_sunM(1:3); 
    v_sunM = X_sunM(4:6); 

    % get sun-Earth state 
    target = 'Earth';
    X_sunE = spice_state(et, target, frame, abcorr, observer); 
    
    % save 
    X_sunM_hist = [X_sunM_hist; X_sunM]; 
    X_sunE_hist = [X_sunE_hist; X_sunE]; 
    et_hist = [et_hist; et]; 
    

    % get angle 
    phi_r = acosd( dot(r_sunE, r_sunM) / (norm(r_sunE)*norm(r_sunM)) ); 
    phi_v = acosd( dot(v_sunE, v_sunM) / (norm(v_sunE)*norm(v_sunM)) ); 
    
end 

X_dep = X_sunE_hist(1,:); 
X_arr = X_sunM_hist(end,:); 

%% simplest case possible 

% Arrival (Mars), AU units 
r_norm_i = norm(X_dep(1:3)); 
r_norm_f  = norm(X_arr(1:3)); 

% velocities 
v_norm_i = sqrt( mu/r_norm_i );
v_norm_f  = sqrt( mu/r_norm_f );

% initial vector 
r_dep = [1 0 0] * r_norm_i; 
v_dep = [0 1 0] * v_norm_i; 
X_dep = [r_dep v_dep]; 

% final vector 
r_arr = [-1 0 0] * r_norm_f; 
v_arr = [0 -1 0] * v_norm_f;
X_arr = [r_arr v_arr]; 


%% HOHMANN TRANSFER 

% Arrival (Mars), AU units 
r_norm_i = norm(X_dep(1:3)); 
r_norm_f  = norm(X_arr(1:3)); 

a_trans = (r_norm_f + r_norm_i)/2; 

v_init = sqrt( mu/r_norm_i );
v_fin  = sqrt( mu/r_norm_f );

% delta v magnitude 
v_trans_a = sqrt( 2*mu/r_norm_i - mu/a_trans ); 
v_trans_b = sqrt( 2*mu/r_norm_f - mu/a_trans ); 

dv_a = v_trans_a - v_init; 
dv_b = v_fin - v_trans_b; 
dv   = norm(dv_a) + norm(dv_b); 

% transfer time 
tau_trans = pi * sqrt( a_trans^3 / mu ); 

% delta v direction 
dv_init = X_dep(4:6) / norm(X_dep(4:6)) * v_trans_a; 
dv_fin  = X_arr(4:6) / norm(X_arr(4:6)) * v_trans_b; 

% initial satellite state for Hohmann transfer 
rv0_sat = X_dep; 
rv0_sat(4:6) = dv_init; 

%% propagate 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% Set run state 
disp('Running sim ...')

% propagate orbit 
% dt = 86400/400; 
t = linspace(et_t0, et_t0 + tau_trans, length(et_hist)); 
tic
[et, X_sunsat] = ode45(@fn.EOM, t, rv0_sat, options); 
toc 

% GMs 
mu_E = 3.986004418e5; 
mu_sun = mu; 
mu_M = 0.042828e6; 

% initialize 
a_Esat_hist = []; 
a_Msat_hist = []; 
a_sunsat_hist = []; 
a_pert_Esun_hist = []; 
a_pert_Msun_hist = []; 

const_vec = 0; 

% determine gravity 
for i = 1:length(X_sunsat) 
    
    % current sun-sat vector 
    r_sunsat = X_sunsat(i, 1:3);
    r_satsun = -r_sunsat; 
    
    if const_vec == 1
        r_sunE = X_dep(1:3);   % constant vector 
        r_sunM = X_arr(1:3); 
    else
        % Current sun-Earth vector 
        target  = 'Earth';
        X_sunE  = spice_state(et(i), target, frame, abcorr, observer); 
        r_sunE = X_sunE(1:3); 
        % Current sun-Mars vector 
        target = 'Mars'; 
        X_sunM = spice_state(et(i), target, frame, abcorr, observer); 
        r_sunM = X_sunM(1:3); 
    end 
        
    % Earth to satellite vector 
    r_Esun = -r_sunE;    
    r_Esat = r_Esun + r_sunsat; 

    % central body accel Earth-satellite
    a_Esat = - mu_E * r_Esat / norm(r_Esat)^3; 
    a_Esat_hist = [a_Esat_hist; norm(a_Esat)]; 
    
    % central body accel sun-satellite 
    a_sunsat = - mu_sun * r_sunsat / norm(r_sunsat)^3; 
    a_sunsat_hist = [a_sunsat_hist; norm(a_sunsat)]; 
    
    % disturbance (third body) Earth-sun 
%     a_dist = - mu_E * r_Esat / norm(r_Esat)^3 - ... 
%         mu_sun * ( r_satsun/norm(r_satsun)^3 + r_Esun/norm(r_Esun)^3); 
    % Vallado disturbance 
    a_pert_Esun = - mu_sun * ( r_satsun/norm(r_satsun)^3 - r_Esun/norm(r_Esun)^3); 
    a_pert_Esun_hist = [a_pert_Esun_hist; norm(a_pert_Esun)]; 
    
    % Mars to satellite vector 
    r_Msun = -r_sunM; 
    r_Msat = r_Msun + r_sunsat;
    
    % central body accel Mars-satellite 
    a_Msat = - mu_M * r_Msat / norm(r_Msat)^3; 
    a_Msat_hist = [a_Msat_hist; norm(a_Msat)]; 
    
    % disturbance (third body) Mars-sun 
    a_pert_Msun = - mu_sun * ( r_satsun/norm(r_satsun)^3 - r_Msun/norm(r_Msun)^3 );
    a_pert_Msun_hist = [a_pert_Msun_hist; norm(a_pert_Msun)]; 
    
end 

t_days = (et-et_t0)/86400; 
% plot accelerations 
figure
    semilogy(t_days, a_Esat_hist, 'linewidth', 2); 
    hold on; grid on; 
    semilogy(t_days, a_sunsat_hist, '--','linewidth', 2); 
    semilogy(t_days, a_pert_Esun_hist, 'linewidth', 1.2); 
%         semilogy((et-et_t0)/86400, a_dist_Esun_P_hist, '-^'); 
    semilogy(t_days, a_Msat_hist, ':', 'linewidth', 1.2)
    semilogy(t_days, a_pert_Msun_hist, '--', 'linewidth', 1.2); 

    legend('2body E-sat', '2body sun-sat', 'pert E-sun (3rd body)', ... 
        '2body M-sat', 'pert M-sun (3rd body)', 'location', 'best')
    ylabel('km/s^2')
    xlabel('Days') 
    title('Gravitational Acceleration')
        
%% SOI crossings 

% find index of acceleration crossing 
da = abs(a_Esat_hist - a_pert_Esun_hist); 
i_min = find(da == min(da)); 

% satellite SOI crossing Earth-sun 
X_sunsat_ES_SOI = X_sunsat(i_min, :); 
r_sunsat = X_sunsat(i_min, 1:3); 

% obtain state at this crossing 
et_i_min = et(i_min); 

% sun-Earth vector 
target = 'Earth';
X_sunE = spice_state(et_i_min, target, frame, abcorr, observer); 
r_sunE = X_sunE(1:3); 
r_Esun = -r_sunE; 

% SOI crossing Earth-sun Earth-sat 
r_SOI_Esun_calc = r_Esun + r_sunsat; 

% mass 
m_E = 5.9724e24; 
m_sun = 1988500e24; 
m_M = 0.64169e24;

% semimajor axis 
a_E = 149.598e6; 
a_M = 227.956e6; 

r_SOI_Esun_ana = ( m_E/m_sun )^(2/5)*a_E; 
r_SOI_Msun_ana = ( m_M/m_sun )^(2/5)*a_M; 

dt = et(2) - et(1); 

for i = 1:length(a_Msat_hist)
    
    % Earth central-perturbation body 
    ratio_Esun(i,:) = a_Esat_hist(i) / a_pert_Esun_hist(i); 
    
    % Mars central-perturbation body 
    ratio_Msun(i,:) = a_Msat_hist(i) / a_pert_Msun_hist(i); 
    
    if i > 1
        dratio_Esun(i,:) = (ratio_Esun(i,:) - ratio_Esun(i-1,:))/dt; 
        dratio_Msun(i,:) = (ratio_Msun(i,:) - ratio_Msun(i-1,:))/dt; 
    else
        dratio_Esun(i,:) = 0; 
        dratio_Msun(i,:) = 0; 
    end 
    
end 

% plot rate of ratio change 
figure
    subplot(2,1,1) 
        semilogy(t_days, dratio_Esun, 'linewidth', 2); 
        grid on; 
        title('Earth to sun ratio') 
    subplot(2,1,2) 
        semilogy(t_days, dratio_Msun, '--', 'linewidth', 2); 
        grid on; 
        title('Mars to sun ratio')
        


%% plot 

leg_hist = []; 
figure()

    % departure 
    leg = quiver3(0, 0, 0, X_dep(1), X_dep(2), X_dep(3)); hold on; grid on; 
    leg_hist = [leg_hist; leg]; 
    
    % arrival 
    leg = quiver3(0, 0, 0, X_arr(1), X_arr(2), X_arr(3)); 
    leg_hist = [leg_hist; leg]; 
    
    % Mars traj 
    leg = plot3(X_sunM_hist(:,1), X_sunM_hist(:,2), X_sunM_hist(:,3), 'r--'); 
    leg_hist = [leg_hist; leg]; 
    plot3(X_sunM_hist(end,1), X_sunM_hist(end,2), X_sunM_hist(end,3), 'r^'); 

    % Earth traj 
    leg = plot3(X_sunE_hist(:,1), X_sunE_hist(:,2), X_sunE_hist(:,3), 'b--'); 
    leg_hist = [leg_hist; leg]; 
    plot3(X_sunE_hist(end,1), X_sunE_hist(end,2), X_sunE_hist(end,3), 'b^'); 
    
    % sun 
    leg = scatter3(0, 0, 0, 'filled'); 
    leg_hist = [leg_hist; leg]; 

    % Hohmann 
    leg = plot3(X_sunsat(:,1), X_sunsat(:,2), X_sunsat(:,3), 'g', 'linewidth', 2); 
    leg_hist = [leg_hist; leg]; 
    plot3(X_sunsat(1,1), X_sunsat(1,2), X_sunsat(1,3), 'go', 'linewidth', 2); 
    plot3(X_sunsat(end,1), X_sunsat(end,2), X_sunsat(end,3), 'g^', 'linewidth', 2); 
    
    % Earth-sun SOI crossing 
    quiver3(0, 0, 0, X_sunsat_ES_SOI(1), X_sunsat_ES_SOI(2), X_sunsat_ES_SOI(3)); 
    
    legend(leg_hist, 'Earth', 'Mars', 'Mars traj', 'Earth traj', 'sun', 'Hohmann' )
    
    axis equal 








