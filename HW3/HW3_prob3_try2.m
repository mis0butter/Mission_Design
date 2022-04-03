%% HW 3
% Junette Hsin 

% sun mu (m^3/s^2)
mu_sun_m = 1.32712440018e20; 
mu_sun_km = mu_sun_m / (1000^3); 
mu = mu_sun_km; 

pos = [100 100 700 700]; 

%% simplest case possible 

% % Arrival (Mars), AU units 
% r_norm_i = norm(X_sunE_hist(1:3)); 
% r_norm_f  = norm(X_sunM_hist(1:3)); 

% semimajor axis 
a_E = 1.4959787e11/1000; 
a_M = 227.956e6; 
r_norm_i = a_E; 
r_norm_f = a_M; 

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

% delta time 
dt = tau_trans / 20000; 
% dt = 100; 

% propagate satellite orbit 
[t, X_sunsat_hist] = ode45(@fn.EOM, [0 : dt : tau_trans], rv0_sat, options); 

% back-propagate Mars from arrival 
[t, X_sunM_hist] = ode45(@fn.EOM, [0 : -dt : -tau_trans], X_arr, options); 
X_sunM_hist = flip(X_sunM_hist); 

% forward propagate Earth 
[t, X_sunE_hist] = ode45(@fn.EOM, [0 : dt : tau_trans], X_dep, options); 

% a useful vector 
X_Esat_hist = -X_sunE_hist + X_sunsat_hist; 
X_satE_hist = -X_Esat_hist; 
X_Msat_hist = -X_sunM_hist + X_sunsat_hist; 
X_satM_hist = -X_Msat_hist; 

test_plot = 0; 
% test if sat to Earth/Mars is correct 
if test_plot == 1
    
    figure()
    
    plot3_xyz(X_sunsat_hist, 'g', 1.2)
    hold on; grid on; 
    plot3_xyz(X_sunsat_hist + X_satE_hist, 'b--', 1.2); 
    plot3_xyz(X_sunsat_hist + X_satM_hist, 'r--', 1.2); 

    plot3( [X_sunsat_hist(1,1), X_sunsat_hist(1,1) + X_satM_hist(1,1)], ... 
        [X_sunsat_hist(1,2), X_sunsat_hist(1,2) + X_satM_hist(1,2)], ... 
        [X_sunsat_hist(1,3), X_sunsat_hist(1,3) + X_satM_hist(1,3)] ); 
    
    plot3( [X_sunsat_hist(end,1), X_sunsat_hist(end,1) + X_satE_hist(end,1)], ... 
        [X_sunsat_hist(end,2), X_sunsat_hist(end,2) + X_satE_hist(end,2)], ... 
        [X_sunsat_hist(end,3), X_sunsat_hist(end,3) + X_satE_hist(end,3)] ); 
    
end 

%% Hohmann plot 

leg_hist = []; 
ftitle = 'Problem 3 - Hohmann Transfer'; 
figure('name', ftitle, 'position', [50 50 700 500])

    % departure 
    leg = quiver3(0, 0, 0, X_dep(1), X_dep(2), X_dep(3)); hold on; grid on; 
    leg_hist = [leg_hist; leg]; 
    
    % arrival 
    leg = quiver3(0, 0, 0, X_arr(1), X_arr(2), X_arr(3)); 
    leg_hist = [leg_hist; leg]; 
    
    % Mars traj 
    leg = plot3_xyz(X_sunM_hist, 'r', 1); 
    leg_hist = [leg_hist; leg]; 
    plot3_xyz(X_sunM_hist, 'ro', 1, 1); 
    plot3_xyz(X_sunM_hist, 'r^', 1, 'end'); 

    % Earth traj 
    leg = plot3_xyz(X_sunE_hist, 'b', 1); 
    leg_hist = [leg_hist; leg]; 
    plot3_xyz(X_sunE_hist, 'bo', 1, 1); 
    plot3_xyz(X_sunE_hist, 'b^', 1, 'end'); 
    
    % sun 
    leg = scatter3(0, 0, 0, 'filled'); 
    leg_hist = [leg_hist; leg]; 

    % Hohmann 
    leg = plot3_xyz(X_sunsat_hist, 'g', 2); 
    leg_hist = [leg_hist; leg]; 
    plot3_xyz(X_sunsat_hist, 'go', 2, 1); 
    plot3_xyz(X_sunsat_hist, 'g^', 2, 'end'); 
    
    legend(leg_hist, 'Earth_{init}', 'Mars_{fin}', 'Mars traj', 'Earth traj', 'sun', 'Hohmann' )
    xlabel('x (km)')
    ylabel('y (km)') 
    zlabel('z (km)') 
    
    view(0, 90)
    
    title(ftitle)
    
    axis equal 


%% gravity 

% GMs 
mu_E = 3.986004418e5; 
mu_sun = mu; 
mu_M = 0.042828e6; 

% mass 
m_E = 5.9724e24; 
m_sun = 1988500e24; 
m_M = 0.64169e24;

% SOI 
r_SOI_Esun_ana = ( m_E/m_sun )^(2/5)*a_E; 
r_SOI_Msun_ana = ( m_M/m_sun )^(2/5)*a_M; 
disp('Analytical r_SOI_Esun norm: ') 
disp(norm(r_SOI_Esun_ana)) 
disp('Analytical r_SOI_Msun norm: ') 
disp(norm(r_SOI_Msun_ana)) 

% initialize 
a_Esat_hist = []; 
a_Msat_hist = []; 
a_sunsat_hist = []; 
a_pert_Esun_hist = []; 
a_pert_sunE_hist = []; 
a_pert_Msun_hist = []; 
a_pert_sunM_hist = []; 

const_vec = 0; 

% determine gravity 
for i = 1:length(X_sunsat_hist) 
    
    % Get current states and positions 
    X_sunsat = X_sunsat_hist(i,:); 
    X_sunE   = X_sunE_hist(i,:); 
    X_sunM   = X_sunM_hist(i,:); 

    r_sunsat = X_sunsat(1:3); 
    r_satsun = -r_sunsat; 
    r_sunE = X_sunE(1:3); 
    r_sunM = X_sunM(1:3); 
    
    % Earth to satellite vector 
    r_Esun = -r_sunE;    
    r_Esat = r_Esun + r_sunsat; 
    r_satE = -r_Esat; 
    
    % Mars to satellite vector 
    r_Msun = -r_sunM; 
    r_Msat = r_Msun + r_sunsat; 
    r_satM = -r_Msat; 

    % central body accel Earth-satellite
    a_Esat = - mu_E * r_Esat / norm(r_Esat)^3; 
    a_Esat_hist = [a_Esat_hist; norm(a_Esat)]; 
    
    % central body accel sun-satellite 
    a_sunsat = - mu_sun * r_sunsat / norm(r_sunsat)^3; 
    a_sunsat_hist = [a_sunsat_hist; norm(a_sunsat)]; 
    
    % disturbance (third body) Earth-sun 
%     a_pert = - mu_E * r_Esat / norm(r_Esat)^3 - ... 
%         mu_sun * ( r_satsun/norm(r_satsun)^3 + r_Esun/norm(r_Esun)^3); 
    % Vallado disturbance 
    % 3rd body pert: Earth-centered, sun pert 
    a_pert_Esun = - mu_sun * ( r_satsun/norm(r_satsun)^3 - r_Esun/norm(r_Esun)^3); 
    a_pert_Esun_hist = [a_pert_Esun_hist; norm(a_pert_Esun)]; 

    % 3rd body pert: sun-centered, Earth pert 
    a_pert_sunE = - mu_E * ( r_satE/norm(r_satE)^3 - r_sunE/norm(r_sunE)^3 ); 
    a_pert_sunE_hist = [a_pert_sunE_hist; norm(a_pert_sunE)]; 
    
    % Mars to satellite vector 
    r_Msun = -r_sunM; 
    r_Msat = r_Msun + r_sunsat;
    
    % central body accel Mars-satellite 
    a_Msat = - mu_M * r_Msat / norm(r_Msat)^3; 
    a_Msat_hist = [a_Msat_hist; norm(a_Msat)]; 
    
    % 3rd-body pert: Mars-centered, sun pert 
    a_pert_Msun = - mu_sun * ( r_satsun/norm(r_satsun)^3 - r_Msun/norm(r_Msun)^3 );
    a_pert_Msun_hist = [a_pert_Msun_hist; norm(a_pert_Msun)]; 
    
    % 3rd-body pert: sun-centered, Mars pert 
    a_pert_sunM = - mu_M * ( r_satM/norm(r_satM)^3 - r_sunM/norm(r_sunM)^3 ); 
    a_pert_sunM_hist = [a_pert_sunM_hist; norm(a_pert_sunM)]; 

    
end 
        
        
%% SOI crossings 

dt = t(2) - t(1); 
for i = 1:length(a_Msat_hist)
    
    % Earth central-perturbation body 
    ratio_Esun(i,:) = a_Esat_hist(i) / a_pert_Esun_hist(i); 
    ratio_sunE(i,:) = a_sunsat_hist(i) / a_pert_sunE_hist(i); 
    
    % Earth-sat norm 
    r_Esat_norm(i,:) = norm(X_Esat_hist(i, 1:3)); 
    
    % Mars central-perturbation body 
    ratio_Msun(i,:) = a_Msat_hist(i) / a_pert_Msun_hist(i); 
    ratio_sunM(i,:) = a_sunsat_hist(i) / a_pert_sunM_hist(i); 
    
    % Mars-sat norm 
    r_Msat_norm(i,:) = norm(X_Msat_hist(i, 1:3)); 
    
    if i > 1
        dratio_Esun(i,:) = norm(ratio_Esun(i,:) - ratio_Esun(i-1,:))/dt; 
        dratio_sunE(i,:) = norm(ratio_sunE(i,:) - ratio_sunE(i-1,:))/dt; 
        dratio_Msun(i,:) = norm(ratio_Msun(i,:) - ratio_Msun(i-1,:))/dt; 
        dratio_sunM(i,:) = norm(ratio_sunM(i,:) - ratio_sunM(i-1,:))/dt; 
    else
        dratio_Esun(i,:) = 0; 
        dratio_sunE(i,:) = 0; 
        dratio_Msun(i,:) = 0; 
        dratio_sunM(i,:) = 0; 
    end 
    
end 

dratio_Esun = abs(ratio_Esun - ratio_sunE); 
i_min = find(dratio_Esun == min(dratio_Esun)); 
t_i_min_Esun = t(i_min); 

dratio_Msun = abs(ratio_Msun - ratio_sunM); 
i_min = find(dratio_Msun == min(dratio_Msun)); 
t_i_min_Msun = t(i_min); 


t_days = t/86400; 

ftitle = 'Problem 3 - Gravitational Acceleration'; 
% plot accelerations 
pos = pos + [50 0 0 0]; 
figure('name', ftitle, 'position', pos)

    subplot(2,1,1) 
    
        semilogy(t_days, a_Esat_hist, 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, a_sunsat_hist, '--','linewidth', 1.2); 
        semilogy(t_days, a_pert_Esun_hist, 'linewidth', 1.2); 
        semilogy(t_days, a_pert_sunE_hist, '--', 'linewidth', 1.2); 
%         semilogy((et-et_t0)/86400, a_dist_Esun_P_hist, '-^'); 
        xline(t_i_min_Esun / 86400, 'b--', 'linewidth', 1.2); 
        
        legend('2body E-sat', '2body sun-sat', 'pert: E central, sun 3rd', ... 
            'pert: sun central, E 3rd', 'Earth SOI', 'location', 'best')
        ylabel('km/s^2')
        title('Earth Gravitational Acceleration')

%         % Earth SOI 
%         xlim([0, t_i_min_Esun/86400 + 1])
        
    subplot(2,1,2) 
        semilogy(t_days, a_Msat_hist, 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, a_sunsat_hist, '--','linewidth', 1.2); 
        semilogy(t_days, a_pert_Msun_hist, 'linewidth', 1.2); 
        semilogy(t_days, a_pert_sunM_hist, '--', 'linewidth', 1.2); 
        xline(t_i_min_Msun / 86400, 'b--', 'linewidth', 1.2); 
        
        legend('2body M-sat', '2body sun-sat', 'pert: M central, 3rd ', ... 
            'pert: sun central, M 3rd', 'Mars SOI', 'location', 'best')
        ylabel('km/s^2')
        title('Mars Gravitational Acceleration'); 
        
%         % Mars SOI 
%         xlim([t_i_min_Msun/86400 - 1, t_days(end) + 1])
        
        sgtitle(ftitle); 
    
        xlabel('Days') 
        

ftitle = 'Problem 3 - Gravitational Acceleration - SOI Vicinity'; 
% plot accelerations 
pos = pos + [50 0 0 0]; 
figure('name', ftitle, 'position', pos)

    subplot(2,1,1) 
    
        semilogy(t_days, a_Esat_hist, 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, a_sunsat_hist, '--','linewidth', 1.2); 
        semilogy(t_days, a_pert_Esun_hist, 'linewidth', 1.2); 
        semilogy(t_days, a_pert_sunE_hist, '--', 'linewidth', 1.2); 
%         semilogy((et-et_t0)/86400, a_dist_Esun_P_hist, '-^'); 
        xline(t_i_min_Esun / 86400, 'b--', 'linewidth', 1.2); 
        
        legend('2body E-sat', '2body sun-sat', 'pert: E central, sun 3rd', ... 
            'pert: sun central, E 3rd', 'Earth SOI', 'location', 'best')
        ylabel('km/s^2')
        title('Earth Gravitational Acceleration')

        % Earth SOI 
        xlim([0, t_i_min_Esun/86400 + 1])
        
    subplot(2,1,2) 
        semilogy(t_days, a_Msat_hist, 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, a_sunsat_hist, '--','linewidth', 1.2); 
        semilogy(t_days, a_pert_Msun_hist, 'linewidth', 1.2); 
        semilogy(t_days, a_pert_sunM_hist, '--', 'linewidth', 1.2); 
        xline(t_i_min_Msun / 86400, 'b--', 'linewidth', 1.2); 
        
        legend('2body M-sat', '2body sun-sat', 'pert: M central, 3rd ', ... 
            'pert: sun central, M 3rd', 'Mars SOI', 'location', 'best')
        ylabel('km/s^2')
        title('Mars Gravitational Acceleration'); 
        
        % Mars SOI 
        xlim([t_i_min_Msun/86400 - 1, t_days(end) + 1])
        
        sgtitle(ftitle); 
    
        xlabel('Days') 

% ------------------------------------------------------------------------
ftitle = 'Problem 3 - Earth, Mars, and Sun acceleration ratios'; 
% plot rate of ratio change 
pos = pos + [50 0 0 0]; 
figure('name', ftitle, 'position', pos)

    % EARTH-SUN RATIO 
    subplot(2,2,1) 
        semilogy(t_days, ratio_Esun, 'b', 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, ratio_sunE, 'k', 'linewidth', 1.2); 

        xline(t_i_min_Esun / 86400, 'b--', 'linewidth', 1.2); 
        
        ylim('auto') 
        legend('E central', 'Sun central', 'E SOI')
        title('Earth to sun ratio') 

%         % Earth SOI 
%         xlim([0, t_i_min_Esun/86400 + 1])

    % MARS-SUN RATIO 
    subplot(2,2,2) 
        semilogy(t_days, ratio_Msun, 'r', 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, ratio_sunM, 'k', 'linewidth', 1.2); 

        xline(t_i_min_Msun / 86400, 'r--', 'linewidth', 1.2); 

        legend('M central', 'Sun central', 'M SOI')
        title('Mars to sun ratio')
        
%         % Mars SOI 
%         xlim([t_i_min_Msun/86400 - 1, t_days(end) + 1])
        
    % EARTH-SUN RATIO CHANGE 
    subplot(2,2,3) 
        semilogy(t_days, dratio_Esun, 'b', 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, dratio_sunE, 'k', 'linewidth', 1.2); 
        ylim('auto') 
        
        xline(t_i_min_Esun / 86400, 'b--', 'linewidth', 1.2); 
        
        legend('E central', 'Sun central', 'E SOI')
        title('Earth to sun ratio rate change magnitude') 
    xlabel('Time (days)') 

%         % Earth SOI 
%         xlim([0, t_i_min_Esun/86400 + 1])

    % MARS-SUN RATIO CHANGE 
    subplot(2,2,4) 
        semilogy(t_days, dratio_Msun, 'r', 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, dratio_sunM, 'k', 'linewidth', 1.2); 
        xline(t_i_min_Msun / 86400, 'r--', 'linewidth', 1.2); 
        
        legend('M central', 'Sun central', 'M SOI')
        title('Mars to sun ratio rate change magnitude') 


%         % Mars SOI 
%         xlim([t_i_min_Msun/86400 - 1, t_days(end) + 1])
        
    sgtitle(ftitle); 
        
    xlabel('Time (days)') 

% ------------------------------------------------------------------------
ftitle = {'Problem 3 - Earth, Mars, and Sun acceleration ratios - SOI vicinity'}; 
% plot rate of ratio change 
pos = pos + [50 0 0 0]; 
figure('name', ftitle{1}, 'position', pos)

    % EARTH-SUN RATIO 
    subplot(2,2,1) 
        semilogy(t_days, ratio_Esun, 'b', 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, ratio_sunE, 'k', 'linewidth', 1.2); 

        xline(t_i_min_Esun / 86400, 'b--', 'linewidth', 1.2); 
        
        ylim('auto') 
        legend('E central', 'Sun central', 'E SOI')
        title('Earth to sun ratio') 

        % Earth SOI 
        xlim([0, t_i_min_Esun/86400 + 1])

    % MARS-SUN RATIO 
    subplot(2,2,2) 
        semilogy(t_days, ratio_Msun, 'r', 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, ratio_sunM, 'k', 'linewidth', 1.2); 

        xline(t_i_min_Msun / 86400, 'r--', 'linewidth', 1.2); 

        legend('M central', 'Sun central', 'M SOI')
        title('Mars to sun ratio')
        
        % Mars SOI 
        xlim([t_i_min_Msun/86400 - 1, t_days(end) + 1])
        
    % EARTH-SUN RATIO CHANGE 
    subplot(2,2,3) 
        semilogy(t_days, dratio_Esun, 'b', 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, dratio_sunE, 'k', 'linewidth', 1.2); 
        ylim('auto') 
        
        xline(t_i_min_Esun / 86400, 'b--', 'linewidth', 1.2); 
        
        legend('E central', 'Sun central', 'E SOI')
        title('Earth to sun ratio rate change magnitude') 
    xlabel('Time (days)') 

        % Earth SOI 
        xlim([0, t_i_min_Esun/86400 + 1])

    % MARS-SUN RATIO CHANGE 
    subplot(2,2,4) 
        semilogy(t_days, dratio_Msun, 'r', 'linewidth', 1.2); 
        hold on; grid on; 
        semilogy(t_days, dratio_sunM, 'k', 'linewidth', 1.2); 
        xline(t_i_min_Msun / 86400, 'r--', 'linewidth', 1.2); 
        
        legend('M central', 'Sun central', 'M SOI')
        title('Mars to sun ratio rate change magnitude') 


        % Mars SOI 
        xlim([t_i_min_Msun/86400 - 1, t_days(end) + 1])
        
    sgtitle(ftitle); 
        
    xlabel('Time (days)') 

%% flight path angle 

for i = 1:length(X_sunsat_hist)
    
    % FLIGHT PATH ANGLE - SUN 
    phi_fpa_sunsat(i,:) = find_fpa(X_sunsat_hist, i); 
        
    % FLIGHT PATH ANGLE - EARTH 
    phi_fpa_Esat(i,:) = find_fpa(X_Esat_hist, i); 
        
    % FLIGHT PATH ANGLE - MARS
    phi_fpa_Msat(i,:) = find_fpa(X_Msat_hist, i); 
    
    % test stuff 
    phi_fpa_sunM(i,:) = find_fpa(X_sunM_hist, i); 
    phi_fpa_sunE(i,:) = find_fpa(X_sunE_hist, i); 
    
end 

ftitle = 'Problem 3 - Flight Path Angle'; 
pos = pos + [50 0 0 0]; 
figure('name', ftitle, 'position', pos)

    subplot(3,1,1) 
        plot(t_days, phi_fpa_sunsat, 'k', 'linewidth', 1.2) 
        hold on; grid on; 
        plot(t_days, phi_fpa_Esat, 'b', 'linewidth', 1.2) 
        plot(t_days, phi_fpa_Msat, 'r', 'linewidth', 1.2) 
        xline(t_i_min_Esun / 86400, 'b--', 'linewidth', 1.2); 
        xline(t_i_min_Msun / 86400, 'r--', 'linewidth', 1.2); 
  
        legend('\phi_{fpa, sun}', '\phi_{fpa, Earth}', '\phi_{fpa, Mars}', ... 
            'Earth SOI', 'Mars SOI')

        xlabel('Time (days)') 
        ylabel('deg') 
        title(ftitle) 

    subplot(3,1,2) 
        plot(t_days, phi_fpa_sunsat, 'k', 'linewidth', 1.2) 
        hold on; grid on; 
        plot(t_days, phi_fpa_Esat, 'b', 'linewidth', 1.2) 
        plot(t_days, phi_fpa_Msat, 'r', 'linewidth', 1.2) 
        xline(t_i_min_Esun / 86400, 'b--', 'linewidth', 1.2); 
        xline(t_i_min_Msun / 86400, 'r--', 'linewidth', 1.2); 

        % Earth SOI 
        xlim([0, t_i_min_Esun/86400 + 1])

        xlabel('Time (days)') 
        ylabel('deg') 
        title(ftitle) 

    subplot(3,1,3) 
        plot(t_days, phi_fpa_sunsat, 'k', 'linewidth', 1.2) 
        hold on; grid on; 
        plot(t_days, phi_fpa_Esat, 'b', 'linewidth', 1.2) 
        plot(t_days, phi_fpa_Msat, 'r', 'linewidth', 1.2) 
        xline(t_i_min_Esun / 86400, 'b--', 'linewidth', 1.2); 
        xline(t_i_min_Msun / 86400, 'r--', 'linewidth', 1.2); 

        % Mars SOI 
        xlim([t_i_min_Msun/86400 - 1, t_days(end) + 1])
        
        fn.move_legend 
        xlabel('Time (days)') 
        ylabel('deg') 
        title(ftitle) 
        
function phi_fpa = find_fpa(X_bodsat_hist, i)
    
    % position unit vector 
        r_sunsat = X_bodsat_hist(i, 1:3); 
        r_sunsat = r_sunsat / norm(r_sunsat); 
    % velocity unit vector 
        v_sunsat = X_bodsat_hist(i, 4:6); 
        v_sunsat = v_sunsat / norm(v_sunsat); 
    % orbit normal 
        h = cross(r_sunsat, v_sunsat); 
    % calculate transverse direction 
        v_transv = cross(h, r_sunsat); 
    % velocity unit vector 
        v_sunsat = X_bodsat_hist(i, 4:6); 
        v_sunsat = v_sunsat / norm(v_sunsat); 
    % fpa 
        phi_fpa = acosd(dot(v_transv, v_sunsat)); 

end 

    
function h = plot3_xyz(X, style, linew, i)

if ~exist('style', 'var') 
    style = ''; 
end 

if ~exist('linew', 'var')
    linew = 1; 
end 

if ~exist('i', 'var')
    h = plot3(X(:,1), X(:,2), X(:,3), style, 'linewidth', linew); 
elseif strcmp(i, 'end')
    h = plot3(X(end,1), X(end,2), X(end,3), style, 'linewidth', linew);     
else
    h = plot3(X(i,1), X(i,2), X(i,3), style, 'linewidth', linew); 
end 

end 


