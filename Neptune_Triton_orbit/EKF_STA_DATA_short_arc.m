% HW 5 
% Junette Hsin 
% TEST EKF SCRIPT 

load WS_ALL_3.mat 

% weighting matrices m --> km, mm --> km 
global R_KJL R_DGO R_ACB 
R_KJL = [(10e-3)^2 0; 0 (0.5e-6)^2]; 
R_DGO = [(5e-3)^2  0; 0 (1e-6)^2]; 
R_ACB = [(1000e-3)^2 0; 0 (0.5e-6)^2]; 

% if DAYS == 1        % 1 day 
%     load LEO_DATA_Apparent.mat 
%     hrs = 24 * 1; 
% elseif DAYS == 3    % 3 days 
%     load('LEO_DATA_Apparent_3Days.mat')    
%     hrs = 24 * 7; 
% else                % 7 days 
% load LEO_DATA_Apparent_3Days.mat 
temp = load('LEO_DATA_Apparent_Days4-6.mat'); 
days4to6 = temp.LEO_DATA_Apparent; 
days4to6(:,2) = days4to6(:,2) + 3*86400; 
% LEO_DATA_Apparent = [LEO_DATA_Apparent; days4to6];   
LEO_DATA_Apparent = days4to6; 
% propagate to last day, which is now the 4th day. hrs = 24*4 = 96 
hrs = 24*4; 
% end

if STATIONS == 0    % use all station data 
elseif STATIONS == 1 % Kwajalein 
    R_DGO = R_DGO * 1e10; 
    R_ACB = R_ACB * 1e10; 
elseif STATIONS == 2 % Diego 
    R_KJL = R_KJL * 1e10; 
    R_ACB = R_ACB * 1e10; 
elseif STATIONS == 3 % STATIONS == 3 % Arecibo
    R_KJL = R_KJL * 1e10; 
    R_DGO = R_DGO * 1e10; 
elseif STATIONS == 4 % Kwajalein and Diego - NO ARECIBO 
    R_ACB = R_ACB * 1e10; 
elseif STATIONS == 5 % Kwajalein and Arecibo - NO DIEGO 
    R_DGO = R_DGO * 1e10; 
else % STATIONS == 6. Diego and Arecibo - NO KWAJALEIN 
    R_KJL = R_KJL * 1e10; 
end

% short-arc overwrites going on 
Yobs_STA   = LEO_DATA_Apparent;
et_obs     = Yobs_STA(:,2) + et_t0; 
% XSTM_prev  = XSTM0_batch; 
X_EKF0     = X_EKF(end,:)'; 
XSTM_prev  = [X_EKF0; STM0];
t_EKF0     = t_X_EKF(end); 
et_t0      = t_EKF0; 
iter       = 0; 
P_prev     = P0; 

% try smaller P0 ... but wider than other cases for short arc 
P_prev     = [ diag([1 1 1]).^2, zeros(3); zeros(3), diag([1e-3 1e-3 1e-3]).^2 ]; 

X_EKF      = [];     t_X_EKF    = [];     Y_prefit   = [];     Y_postfit  = []; 
Lpre_mat   = [];     Lpost_mat  = [];     sigma3_pre = [];     sigma3_post = []; 

% EKF 
tic
for i = 1:length(et_obs) 

    % keep track of iterations 
    iter = iter + 1; 
    sprintf('iter = %d', iter)

    % Propagate state 
    if     i == 1 && et_obs(1) == et_t0; t_prop = et_obs(i); 
    elseif i == 1;                       t_prop = [et_t0 : 60 : et_obs(1) ]; 
    else                                
%         t_prop = [et_obs(i-1) : 60 : et_obs(i)]; 
        t_prop = [t_lt et_obs(i)]; 
    end

    % EKF. All data, range, or range-only 
    [t_XSTM, XSTM, Xstar, Y_pre, Y_post, P, L_pre, L_post, t_lt] = fn.EKF(Yobs_STA, XSTM_prev, nX, ... 
        epochs(1), t_prop, options, Amat_fn, Ht_fn, Ht_bias_fn, Ht_r_fn, Ht_rr_fn, P_prev, DATA); 

    % save states from current iteration 
    if i == 1 && et_obs(1) == et_t0; X_EKF = XSTM(1:nX)'; 
    else;                            X_EKF = [X_EKF; XSTM(:, 1:nX)]; 
    end 

    t_X_EKF   = [t_X_EKF; t_XSTM]; 
    Y_prefit  = [Y_prefit; Y_pre]; 
    Y_postfit = [Y_postfit; Y_post]; 

    % update measurement for next iteration 
    XSTM_prev = [Xstar; STM0]; 
    P_prev    = P;

    % innovations covariance 
    Lpre_mat  = [Lpre_mat; L_pre]; 
    Lpost_mat = [Lpost_mat; L_post]; 

    % 3-sigma STD 
    if DATA == 0
        sigma3_pre  = [sigma3_pre; sqrt(L_pre(1,1))*3, sqrt(L_pre(2,2))*3];  
        sigma3_post = [sigma3_post; sqrt(L_post(1,1))*3, sqrt(L_post(2,2))*3];  
    else  
        sigma3_pre  = [sigma3_pre; sqrt(L_pre(1,1))*3];  
        sigma3_post = [sigma3_post; sqrt(L_post(1,1))*3];  
    end

end 
toc

%% Propagate to last period of time 

t_prop = [ et_obs(end) : 60 : et_t0 + 60*60*hrs + 420.000000357591 ]; 
[t_XSTM, XSTM] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn, nX), t_prop, XSTM_prev, options); 
Xi   = XSTM(end,1:nX)'; 
STMi = XSTM(end,nX+1:end); 
STMi = reshape(STMi, [nX nX]); 

% Time update + process noise 
% dt       = t_prop(end) - t_prop(1); 
dt       = 60; 
Q        = diag( (100e-10)^2 * [1 1 1] ); 
Gamma    = [diag( dt^2/2 * [1 1 1] ); diag([dt dt dt])]; 
P_noise  = Gamma * Q * Gamma'; 
P_bar    = STMi * P_prev * STMi' + P_noise; 

% Save propagated states 
X_prop   = XSTM(:,1:nX); 
t_X_prop = t_XSTM; 

% PLOT RESIDUALS 
EKF_res; 

%% Save data 

save('WS_ALL_7_short_arc'); 


%% radial-intrack-cross-track frame transformation for best estimate 

if plot_RSW == 1
    
    % 1 day best 
    X_best = [-6330.24806665475
               3306.88524247528
               127.818471301279
              -3.43799666089757
              -6.63343176758256
             -0.235596451267133 ]; 

    T_best = fn.ECItoRSW_T(X_best); 

    dr_ECI = X_best(1:3) - Xi(1:3); 

    % transform all measurements 
    dr_RSW = T_best * dr_ECI; 

    % Plot radial-intrack-crosstrack 
    ftitle = 'Radial-Intrack-Crosstrack'; 
    figh = figure('name', ftitle); 
        subplot(3,1,1)
            scatter(dr_RSW(1), dr_RSW(2)); hold on; 
            P = P_bar(1:2, 1:2); 
            h3 = fn.plot_gaussian_ellipsoid([dr_RSW(1) dr_RSW(2)], P); 
            xlabel('R (km)')
            ylabel('S (km)')
            title('Radial-Intrack') 
        subplot(3,1,2) 
            scatter(dr_RSW(1), dr_RSW(3)); hold on; 
            P = [P_bar(1), P_bar(1,3); P_bar(3,1), P_bar(3,3)]; 
            h3 = fn.plot_gaussian_ellipsoid([dr_RSW(1) dr_RSW(3)], P); 
            xlabel('R (km)')
            ylabel('W (km)')
            title('Radial-Crosstrack') 
        subplot(3,1,3) 
            scatter(dr_RSW(2), dr_RSW(3)); hold on; 
            P = P_bar(2:3, 2:3); 
            h3 = fn.plot_gaussian_ellipsoid([dr_RSW(2) dr_RSW(3)], P); 
            xlabel('S (km)')
            ylabel('W (km)')
            title('Intrack-Crosstrack')

            % Vallado ed 4 p. 229 
end


%% Plot satellite position 

if plot_orbit == 1
    ftitle = 'JahSat Orbit'; 
    figure('name', ftitle); 
        plot3(XSTM_ref0(:,1), XSTM_ref0(:,2), XSTM_ref0(:,3)); hold on; grid on; 
        plot3(XSTM_batch(:,1), XSTM_batch(:,2), XSTM_batch(:,3)); 
        plot3(X_EKF(:,1), X_EKF(:,2), X_EKF(:,3));  
        plot3(X_prop(:,1), X_prop(:,2), X_prop(:,3))
        xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)'); 
        legend('initial', 'batch', 'EKF', 'prop') 
        title(ftitle)
end




