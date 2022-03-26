function [t_XSTM, XSTM, X_update, Y_prefit, Y_postfit, P_update, L_pre, L_post, t_lt] = ... 
    EKF(Yobs_STA, XSTM_prev, nX, et_t0, t_prop, options, Amat_fn, Ht_fn, Ht_bias_fn, Ht_r_fn, Ht_rr_fn, P_prev, DATA)

global R_KJL R_DGO R_ACB 
global r_KJL_ECEF r_DGO_ECEF r_ACB_ECEF 

% Integrate ref trajectory and STM from t = i-1 (prev) to t = i (curr) 
if length(t_prop) > 1
    [t_XSTM, XSTM] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn, nX), t_prop, XSTM_prev, options); 
    Xi   = XSTM(end,1:nX)'; 
    STMi = XSTM(end,nX+1:end); 
else
    t_XSTM = t_prop; 
    XSTM = XSTM_prev; 
    Xi   = XSTM(1:nX); 
    STMi = XSTM(nX+1:end); 
end 
STMi = reshape(STMi, [nX nX]); 

% find index for same time state and observation 
t_Y   = Yobs_STA(:,2) + et_t0; % time after initial epoch 
ti_X  = t_XSTM(end); 
i_Y   = find(t_Y == ti_X); 
Yi    = Yobs_STA(i_Y, :); 

% observation covariance 
if Yi(1) == 1 
    R = R_KJL; 
    r_STA_ECEF = r_KJL_ECEF; 
elseif Yi(1) == 2 
    R = R_DGO; 
    r_STA_ECEF = r_DGO_ECEF; 
else 
    R = R_ACB; 
    r_STA_ECEF = r_ACB_ECEF; 
end 

% get JD time 
JD_UTC = cspice_et2utc(t_XSTM(end), 'J', 10); 
JD_UTC = str2num(extractAfter(JD_UTC, 'JD ')); 

% Convert station to ECI frame 
XSi_OG = rv_ECEFtoECI(JD_UTC, r_STA_ECEF, [0; 0; 0]); 
MJD_UTC = JD_UTC - 2400000.5; 
XSi = ECEF2ECI( MJD_UTC, [r_STA_ECEF; 0; 0; 0] ); 

% Time update + process noise 
dt    = t_prop(end) - t_prop(1); 
Q     = diag( (10e-10)^2 * [1 1 1] ); 
Gamma = [diag( dt^2/2 * [1 1 1] ); diag([dt dt dt])]; 
P_noise = Gamma * Q * Gamma'; 
Ppre_bar = STMi * P_prev * STMi' + P_noise; 

% Y prefit 
Y_prefit(1:2) = Yi(1:2); 
Y_prefit(3:4) = fn.G_fn(Xi, XSi)'; 

% PREFIT Observation-state matrix 
if DATA == 0 || DATA == 3
    Hti_pre = Ht_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 
elseif DATA == 1
    R = R(1,1); 
    Hti_pre = Ht_r_fn(Xi(1), Xi(2), Xi(3), XSi(1), XSi(2), XSi(3)); 
else % DATA == 2
    R = R(2,2); 
    Hti_pre = Ht_rr_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 
end

% PREFIT Innovation (information) covariance 
L_pre = (Hti_pre * Ppre_bar * Hti_pre' + R); 

% CORRECT + UPDATE. Light time correction 
c     = 299792.458; % km/s 
lt    = Y_prefit(3) / c;  % range / c = delta time (s) 
t_lt  = ti_X - lt; 
t_back = [ti_X, t_lt];  
if length(t_prop) > 1 
    XSTM_end = XSTM(end,:); 
else
    XSTM_end = XSTM; 
end
nX = length(Xi); 
[~, XSTM_lt] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn, nX), t_back, XSTM_end, options); 
Xi_lt        = XSTM_lt(end,1:nX)'; 

% Satellite ECI coords with aberration and light time corrections. IF
% RANGE-RATE ONLY, TURN OFF ABERRATION 
v_STA_ECI      = XSi(4:6); 
Xi_lt_abr      = Xi_lt; 
% if ~isequal(DATA, 2)
%     Xi_lt_abr(1:3) = Xi_lt(1:3) + lt * v_STA_ECI; 
% end 

if DATA == 0 || DATA == 3

    % Compute Hti again, with satellite ECI lt and abr corrected state, and same station ECI state 
    Hti_lt_abr = Ht_bias_fn(Xi_lt_abr(1), Xi_lt_abr(2), Xi_lt_abr(3), Xi_lt_abr(4), Xi_lt_abr(5), Xi(6), ... 
        XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 
%     Hti_lt_abr = Ht_fn(Xi_lt_abr(1), Xi_lt_abr(2), Xi_lt_abr(3), Xi_lt_abr(4), Xi_lt_abr(5), Xi(6), ... 
%         XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 

    % Gain matrix 
    Ki = Ppre_bar * Hti_lt_abr' * inv( Hti_lt_abr * Ppre_bar * Hti_lt_abr' + R ); 
    
    % Obtain y difference 
    yi = Yi(3:4)' - fn.G_bias_fn(Xi_lt_abr, XSi, Yi(1)); 
%     yi = Yi(3:4)' - fn.G_fn(Xi_lt_abr, XSi); 
%     yi(1) = yi(1) + rbias; 
    
elseif DATA == 1
    
    % Compute Hti again, with satellite ECI lt and abr corrected state, and same station ECI state 
    Hti_lt_abr = Ht_r_fn(Xi_lt_abr(1), Xi_lt_abr(2), Xi_lt_abr(3), XSi(1), XSi(2), XSi(3)); 

    % Gain matrix 
    Ki = Ppre_bar * Hti_lt_abr' * inv( Hti_lt_abr * Ppre_bar * Hti_lt_abr' + R); 

    % Obtain y difference (range only) 
    yi = Yi(3:4)' - fn.G_fn(Xi_lt_abr, XSi); 
%     yi = Yi(3:4)' - fn.G_bias_fn(Xi_lt_abr, XSi, Yi(2)); 
    yi = yi(1); 
    
elseif DATA == 2
    
    % Compute Hti again, with satellite ECI lt and abr corrected state, and same station ECI state 
    Hti_lt_abr = Ht_rr_fn(Xi_lt_abr(1), Xi_lt_abr(2), Xi_lt_abr(3), Xi_lt_abr(4), Xi_lt_abr(5), Xi_lt_abr(6), ... 
        XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 

    % Gain matrix 
    Ki = Ppre_bar * Hti_lt_abr' * inv( Hti_lt_abr * Ppre_bar * Hti_lt_abr' + R); 

    % Obtain y difference (range-rate only) 
    yi = Yi(3:4)' - fn.G_fn(Xi_lt_abr, XSi); 
    yi = yi(2); 

end

% Measurement and reference orbit update. Add xhat to dynamically propagated satellite ECI state  
xhat     = Ki * yi; 
if norm(xhat(1:3)) > 50 
   disp('xhat pos > 50'); return 
end
X_update = Xi_lt + xhat; 
nX       = length(Xi); 
P_update = ( eye(nX) - Ki * Hti_pre ) * Ppre_bar; 

% POSTFIT Observation-state matrix <-- FIX, USE XSi 
if DATA == 0 || DATA == 3
    Hti_post = Ht_fn(X_update(1), X_update(2), X_update(3), X_update(4), X_update(5), X_update(6), ... 
        XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 
elseif DATA == 1
    Hti_post = Ht_r_fn(X_update(1), X_update(2), X_update(3), XSi(1), XSi(2), XSi(3)); 
else
    Hti_post = Ht_rr_fn(X_update(1), X_update(2), X_update(3), X_update(4), X_update(5), X_update(6), ... 
        XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6));     
end

% POSTFIT Innovation (information) covariance 

L_post = (Hti_post * P_update * Hti_post' + R); 

% Y postfit 
Y_postfit(1:2) = Yi(1:2); 
Y_postfit(3:4) = fn.G_fn(X_update, XSi)'; 


end 

%% subfunctions 

function XSi = rv_ECEFtoECI(JD_UTC, r_STA_ECEF, v_STA_ECEF)

global wE 

    r_STA_ECI  = fn.ECEFtoECI(JD_UTC, r_STA_ECEF); 
    v_STA_ECI  = v_STA_ECEF + cross([ 0 0 wE ]', r_STA_ECEF); 
    v_STA_ECI  = fn.ECEFtoECI(JD_UTC, v_STA_ECI); % Technically wrong. Look in Vallado p. 228 
    XSi        = [r_STA_ECI; v_STA_ECI]; 

end

% function [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update, t_lt] = DATA_0_corr(Yi, Xi, XSi, Ht_fn, Ppre_bar, R, ... 
%     ti_X, t_prop, XSTM, options, Amat_fn)
%     
%     % Y prefit 
%     Y_prefit(1:2) = Yi(1:2); 
%     Y_prefit(3:4) = fn.G_fn(Xi, XSi)'; 
% 
%     % PREFIT Observation-state matrix 
%     Hti_pre = Ht_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 
% 
%     % PREFIT Innovation (information) covariance 
%     L_pre = (Hti_pre * Ppre_bar * Hti_pre' + R); 
%     
%     %%%%%% CORRECT + UPDATE 
% 
%     % Light time correction 
%     c     = 299792.458; % km/s 
%     lt    = Y_prefit(3) / c;  % range / c = delta time (s) 
%     t_lt  = ti_X - lt; 
%     t_back = [ti_X, t_lt];  
%     if length(t_prop) > 1 
%         XSTM_end = XSTM(end,:); 
%     else
%         XSTM_end = XSTM; 
%     end
%     nX = length(Xi); 
%     [~, XSTM_lt] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn, nX), t_back, XSTM_end, options); 
%     Xi_lt        = XSTM_lt(end,1:nX)'; 
%     
%     % Satellite ECI coords with aberration and light time corrections 
%     v_STA_ECI      = XSi(4:6); 
%     Xi_lt_abr      = Xi_lt; 
%     Xi_lt_abr(1:3) = Xi_lt(1:3) + lt * v_STA_ECI; 
%     
%     % Compute Hti again, with satellite ECI lt and abr corrected state, and same station ECI state 
%     Hti_lt_abr = Ht_fn(Xi_lt_abr(1), Xi_lt_abr(2), Xi_lt_abr(3), Xi_lt_abr(4), Xi_lt_abr(5), Xi(6), ... 
%         XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 
%     
%     % Gain matrix 
%     Ki = Ppre_bar * Hti_lt_abr' * inv( Hti_lt_abr * Ppre_bar * Hti_lt_abr' + R ); 
%     
%     % Obtain y difference 
% %     yi = Yi(3:4)' - fn.G_bias_fn(Xi_lt_abr, XSi, Yi(2)); 
%     yi = Yi(3:4)' - fn.G_fn(Xi_lt_abr, XSi); 
%     
%     % Measurement and reference orbit update. Add xhat to dynamically propagated satellite ECI state  
%     xhat     = Ki * yi; 
%     X_update = Xi_lt + xhat; 
%     nX       = length(Xi); 
%     P_update = ( eye(nX) - Ki * Hti_pre ) * Ppre_bar; 
%     % Pi    = ( eye(nX) - Ki * Hti ) * Pi_bar * ( eye(nX) - Ki * Hti )' + Ki * R * Ki'; 
% 
%     % POSTFIT Observation-state matrix <-- FIX, USE XSi 
%     Hti_post = Ht_fn(X_update(1), X_update(2), X_update(3), X_update(4), X_update(5), X_update(6), ... 
%         XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 
% 
%     % POSTFIT Innovation (information) covariance 
%     L_post = (Hti_post * P_update * Hti_post' + R); 
% 
%     % Y postfit 
%     Y_postfit(1:2) = Yi(1:2); 
%     Y_postfit(3:4) = fn.G_fn(X_update, XSi)'; 
% 
% end 
% 
% function [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update] = DATA_0(Yi, Xi, XSi, Ht_fn, Pi_bar, R)
% 
%     % Y prefit 
%     Y_prefit(1:2) = Yi(1:2); 
%     Y_prefit(3:4) = fn.G_fn(Xi, XSi)'; 
% 
%     % PREFIT Observation-state matrix 
%     Hti_pre = Ht_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 
% 
%     % PREFIT Innovation (information) covariance 
%     L_pre = (Hti_pre * Pi_bar * Hti_pre' + R); 
% 
%     % Gain matrix 
%     Ki = Pi_bar * Hti_pre' * inv( Hti_pre * Pi_bar * Hti_pre' + R ); 
% 
%     % Obtain y difference 
%     yi = Yi(3:4)' - fn.G_fn(Xi, XSi); 
% 
%     % Measurement and reference orbit update 
%     xhat     = Ki * yi; 
%     X_update = Xi + xhat; 
%     nX       = length(Xi); 
%     P_update = ( eye(nX) - Ki * Hti_pre ) * Pi_bar;
% 
%     % POSTFIT Observation-state matrix 
%     Hti_post = Ht_fn(X_update(1), X_update(2), X_update(3), X_update(4), X_update(5), X_update(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 
% 
%     % Innovation (information) covariance 
%     L_post = (Hti_post * P_update * Hti_post' + R); 
% 
%     % Y postfit 
%     Y_postfit(1:2) = Yi(1:2); 
%     Y_postfit(3:4) = fn.G_fn(X_update, XSi)'; 
% 
% end 
% 
% function [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update] = DATA_1(Yi, Xi, XSi, Ht_r_fn, Pi_bar, R)
% 
%     R = R(1,1); 
% 
%     % Y prefit 
%     Y_prefit(1:2) = Yi(1:2); 
%     y             = fn.G_fn(Xi, XSi)'; 
%     Y_prefit(3)   = y(1); 
% 
%     % PREFIT Observation-state matrix 
%     Hti_pre = Ht_r_fn(Xi(1), Xi(2), Xi(3), XSi(1), XSi(2), XSi(3)); 
% 
%     % PREFIT Innovation (information) covariance 
%     L_pre = (Hti_pre * Pi_bar * Hti_pre' + R); 
% 
%     % Gain matrix 
%     Ki = Pi_bar * Hti_pre' * inv( Hti_pre * Pi_bar * Hti_pre' + R ); 
% 
%     % Obtain y difference (range only) 
%     yi = Yi(3:4)' - fn.G_fn(Xi, XSi); 
%     yi = yi(1); 
% 
%     % Measurement and reference orbit update 
%     xhat     = Ki * yi; 
%     X_update = Xi + xhat; 
%     nX       = length(Xi); 
%     P_update = ( eye(nX) - Ki * Hti_pre ) * Pi_bar;
% 
%     % POSTFIT Observation-state matrix 
%     Hti_post = Ht_r_fn(X_update(1), X_update(2), X_update(3), XSi(1), XSi(2), XSi(3)); 
% 
%     % POSTFIT Innovation (information) covariance 
%     L_post = (Hti_post * P_update * Hti_post' + R); 
% 
%     % Y postfit 
%     Y_postfit(1:2) = Yi(1:2); 
%     y              = fn.G_fn(X_update, XSi)'; 
%     Y_postfit(3)   = y(1); 
% 
% end
% 
% function [X_update, Y_prefit, L_pre, L_post, Y_postfit, P_update] = DATA_2(Yi, Xi, XSi, Ht_rr_fn, Pi_bar, R)
% 
%     R = R(2,2); 
% 
%     % Y prefit 
%     Y_prefit(1:2) = Yi(1:2); 
%     y             = fn.G_fn(Xi, XSi)'; 
%     Y_prefit(3)   = y(2); 
% 
%     % PREFIT Observation-state matrix 
%     Hti_pre = Ht_rr_fn(Xi(1), Xi(2), Xi(3), Xi(4), Xi(5), Xi(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 
% 
%     % PREFIT Innovation (information) covariance 
%     L_pre = (Hti_pre * Pi_bar * Hti_pre' + R); 
% 
%     % Gain matrix 
%     Ki = Pi_bar * Hti_pre' * inv( Hti_pre * Pi_bar * Hti_pre' + R ); 
% 
%     % Obtain y difference (range only) 
%     yi = Yi(3:4)' - fn.G_fn(Xi, XSi); 
%     yi = yi(2); 
% 
%     % Measurement and reference orbit update 
%     xhat     = Ki * yi; 
%     X_update = Xi + xhat; 
%     nX       = length(Xi); 
%     P_update = ( eye(nX) - Ki * Hti_pre ) * Pi_bar;
% 
%     % POSTFIT Observation-state matrix 
%     Hti_post = Ht_rr_fn(X_update(1), X_update(2), X_update(3), X_update(4), X_update(5), X_update(6), XSi(1), XSi(2), XSi(3), XSi(4), XSi(5), XSi(6)); 
% 
%     % Innovation (information) covariance 
%     L_post = (Hti_post * P_update * Hti_post' + R); 
% 
%     % Y postfit 
%     Y_postfit(1:2) = Yi(1:2); 
%     y              = fn.G_fn(X_update, XSi)'; 
%     Y_postfit(3)   = y(2); 
% 
% end
% 
% 
% 


