%% HW 2 
% Junette Hsin 

% close all; 
clear; 
addpath(genpath('Lambert Battin'))

%% transfer angle = 75 deg 

%  Define parameters for a state lookup:
% t0      = 'Oct 20, 2020 11:00 AM CST'; 
t0      = 'May 22, 2020'; 

phi_t_des = 75; 
[ell_1_min, ell_2_min, amin_AU, emin] = lambert_prob(t0, phi_t_des, 1); 


%% phi = 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180 degrees

phi_d_hist  = []; 
phi_a_hist  = []; 
tof_hist = []; 
phi_t_hist = []; 

phi0 = 15; 
for phi = [ phi0 : 15 : 165, 179] 
    
    phi_t_des = phi; 
    [ell_1_min, ell_2_min] = lambert_prob(t0, phi_t_des, 0); 
    
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

figure()
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
        ylabel('deg') 
        legend(h, 'ell 1 short', 'ell 1 long', 'ell 2 short', 'ell 2 long', 'location', 'eastoutside'); 
        
    xlabel('Transfer Angle Phi (deg)') 
    sgtitle('Min Energy Parameters') 

%% save plots 

% savePDF(gcf, get(gcf, 'name'), 'HW2/latex')





