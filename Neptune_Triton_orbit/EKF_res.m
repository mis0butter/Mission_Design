


%% Calculate residuals 

Ypre_KJL   = [];    Ypre_DGO   = [];    Ypre_ACB   = []; 
Ypost_KJL  = [];    Ypost_DGO  = [];    Ypost_ACB  = []; 
sigma3_KJL = [];    sigma3_DGO = [];    sigma3_ACB = []; 

% Extract states that correspond with station measurements 
for i = 1:length(Y_postfit)

    ti  = Y_postfit(i, 2); 
    ti  = ti + et_t0; 
    i_X = find(t_X_EKF == ti); 
    i_STA = Y_postfit(i, 1); 

    if i_STA == 1
        Ypre_KJL   = [Ypre_KJL; Y_prefit(i,:)]; 
        Ypost_KJL  = [Ypost_KJL; Y_postfit(i,:)];  
    elseif i_STA == 2
        Ypre_DGO   = [Ypre_DGO; Y_prefit(i,:)]; 
        Ypost_DGO  = [Ypost_DGO; Y_postfit(i,:)];  
    else
        Ypre_ACB   = [Ypre_ACB; Y_prefit(i,:)]; 
        Ypost_ACB  = [Ypost_ACB; Y_postfit(i,:)];  
    end 

end 


% Station data 
i_STA    = find(LEO_DATA_Apparent(:, 1) == 1); 
Yobs_KJL = LEO_DATA_Apparent(i_STA, :); 
i_STA    = find(LEO_DATA_Apparent(:, 1) == 2); 
Yobs_DGO = LEO_DATA_Apparent(i_STA, :);
i_STA    = find(LEO_DATA_Apparent(:, 1) == 3); 
Yobs_ACB = LEO_DATA_Apparent(i_STA, :);
[dpre_err_KJL, dpre_rms_KJL, vpre_err_KJL, vpre_rms_KJL] = calc_res_all(Yobs_KJL, Ypre_KJL, DATA); 
[dpre_err_DGO, dpre_rms_DGO, vpre_err_DGO, vpre_rms_DGO] = calc_res_all(Yobs_DGO, Ypre_DGO, DATA); 
[dpre_err_ACB, dpre_rms_ACB, vpre_err_ACB, vpre_rms_ACB] = calc_res_all(Yobs_ACB, Ypre_ACB, DATA); 
[dpost_err_KJL, dpost_rms_KJL, vpost_err_KJL, vpost_rms_KJL] = calc_res_all(Yobs_KJL, Ypost_KJL, DATA); 
[dpost_err_DGO, dpost_rms_DGO, vpost_err_DGO, vpost_rms_DGO] = calc_res_all(Yobs_DGO, Ypost_DGO, DATA); 
[dpost_err_ACB, dpost_rms_ACB, vpost_err_ACB, vpost_rms_ACB] = calc_res_all(Yobs_ACB, Ypost_ACB, DATA); 
t_KJL = Yobs_KJL(:,2) / 3600; 
t_DGO = Yobs_DGO(:,2) / 3600; 
t_ACB = Yobs_ACB(:,2) / 3600; 
t_ALL = Y_prefit(:,2) / 3600 ; 

if DATA == 0

    if STATIONS == 0;       ftitle = 'All Stations, All Data Residuals'; 
    elseif STATIONS == 1;   ftitle = 'Kwajalein Only, All Data Residuals'; 
    elseif STATIONS == 2;   ftitle = 'Diego-Garcia Only, All Data Residuals'; 
    elseif STATIONS == 3;   ftitle = 'Arecibo Only, All Data Residuals';    
    elseif STATIONS == 4;   ftitle = 'Kwajalein and Diego-Garcia, NO Arecibo, All Data Residuals'; 
    elseif STATIONS == 5;   ftitle = 'Kwajalein and Arecibo, NO Diego-Garcia, All Data Residuals'; 
    else;                   ftitle = 'Diego-Garcia and Arecibo, NO Kwajalein, All Data Residuals'; % STATIONS == 6  
    end

    % 1 2 3 4 
    % 5 6 7 8 
    % 9 10 11 12 
    % 13 14 15 16 
    figure('name', ftitle); 
        % first row: subplot(4,4,1:4) 
        plot_res(4, 4, 1:4, 'PREFIT range residuals', t_ALL, sigma3_pre(:,1), t_KJL, t_DGO, t_ACB, dpre_err_KJL, ... 
            dpre_err_DGO, dpre_err_ACB, dpre_rms_KJL, dpre_rms_DGO, dpre_rms_ACB, 'km')
        % second row: subplot(4,4,5:8)
        plot_res(4, 4, 5:8, 'PREFIT range-rate residuals', t_ALL, sigma3_pre(:,2), t_KJL, t_DGO, t_ACB, vpre_err_KJL, ... 
            vpre_err_DGO, vpre_err_ACB, vpre_rms_KJL, vpre_rms_DGO, vpre_rms_ACB, 'km/s')
        % third row: subplot(4,4,9:12) 
        plot_res(4, 4, 9:12, 'POSTFIT range residuals', t_ALL, sigma3_post(:,1), t_KJL, t_DGO, t_ACB, dpost_err_KJL, ... 
            dpost_err_DGO, dpost_err_ACB, dpost_rms_KJL, dpost_rms_DGO, dpost_rms_ACB, 'km')
        % fourth row: subplot(4,4,13:16)
        plot_res(4, 4, 13:16, 'POSTFIT range-rate residuals', t_ALL, sigma3_post(:,2), t_KJL, t_DGO, t_ACB, vpost_err_KJL, ... 
            vpost_err_DGO, vpost_err_ACB, vpost_rms_KJL, vpost_rms_DGO, vpost_rms_ACB, 'km/s')

        xlabel('Time (hr)') 
        sgtitle(ftitle); 

elseif DATA == 1

    ftitle = 'All Stations, Range Only Residuals'; 
    figure('name', ftitle); 
    % first row: subplot(4,4,1:4) 
    plot_res(4, 4, 1:4, 'PREFIT range residuals', t_ALL, sigma3_pre, t_KJL, t_DGO, t_ACB, dpre_err_KJL, ... 
        dpre_err_DGO, dpre_err_ACB, dpre_rms_KJL, dpre_rms_DGO, dpre_rms_ACB, 'km')
    % second row: subplot(4,4,5:8)
    plot_res(4, 4, 5:8, 'PREFIT range-rate residuals', t_ALL, NaN * ones(size(t_ALL)), t_KJL, t_DGO, t_ACB, vpre_err_KJL, ... 
        vpre_err_DGO, vpre_err_ACB, vpre_rms_KJL, vpre_rms_DGO, vpre_rms_ACB, 'km/s')
    % third row: subplot(4,4,9:12) 
    plot_res(4, 4, 9:12, 'POSTFIT range residuals', t_ALL, sigma3_post, t_KJL, t_DGO, t_ACB, dpost_err_KJL, ... 
        dpost_err_DGO, dpost_err_ACB, dpost_rms_KJL, dpost_rms_DGO, dpost_rms_ACB, 'km')
    % fourth row: subplot(4,4,13:16)
    plot_res(4, 4, 13:16, 'POSTFIT range-rate residuals', t_ALL, NaN * ones(size(t_ALL)), t_KJL, t_DGO, t_ACB, vpost_err_KJL, ... 
        vpost_err_DGO, vpost_err_ACB, vpost_rms_KJL, vpost_rms_DGO, vpost_rms_ACB, 'km/s')

    xlabel('Time (hr)') 
    sgtitle(ftitle); 
    
else

    ftitle = 'All Stations, Range-Rate Only Residuals'; 
    figure('name', ftitle); 
    % first row: subplot(4,4,1:4) 
    plot_res(4, 4, 1:4, 'PREFIT range residuals', t_ALL, NaN * ones(size(t_ALL)), t_KJL, t_DGO, t_ACB, dpre_err_KJL, ... 
        dpre_err_DGO, dpre_err_ACB, dpre_rms_KJL, dpre_rms_DGO, dpre_rms_ACB, 'km')
    % second row: subplot(4,4,5:8)
    plot_res(4, 4, 5:8, 'PREFIT range-rate residuals', t_ALL, sigma3_pre, t_KJL, t_DGO, t_ACB, vpre_err_KJL, ... 
        vpre_err_DGO, vpre_err_ACB, vpre_rms_KJL, vpre_rms_DGO, vpre_rms_ACB, 'km/s')
    % third row: subplot(4,4,9:12) 
    plot_res(4, 4, 9:12, 'POSTFIT range residuals', t_ALL, NaN * ones(size(t_ALL)), t_KJL, t_DGO, t_ACB, dpost_err_KJL, ... 
        dpost_err_DGO, dpost_err_ACB, dpost_rms_KJL, dpost_rms_DGO, dpost_rms_ACB, 'km')
    % fourth row: subplot(4,4,13:16)
    plot_res(4, 4, 13:16, 'POSTFIT range-rate residuals', t_ALL, sigma3_post, t_KJL, t_DGO, t_ACB, vpost_err_KJL, ... 
        vpost_err_DGO, vpost_err_ACB, vpost_rms_KJL, vpost_rms_DGO, vpost_rms_ACB, 'km/s')

    xlabel('Time (hr)') 
    sgtitle(ftitle); 
    
end



%% Subfunctions 

function plot_res(m, n, ivec, ftitle, t_ALL, sigma3, t_KJL, t_DGO, t_ACB, pre_err_KJL, ... 
    pre_err_DGO, pre_err_ACB, pre_rms_KJL, pre_rms_DGO, pre_rms_ACB, yl)

    if ~isnan(sigma3) 
%         yrange = 3 * mean(real(sigma3)); 
        yrange = 3 * max( [ max(abs( pre_rms_DGO )) max(abs( pre_rms_KJL )) max(abs( pre_rms_ACB )) ] ); 
        
    end
    
    subplot(m, n, ivec(end) ) 
        stext = { sprintf('KJL mean = %.3g', mean(pre_err_KJL)); sprintf('KJL RMS = %.3g', pre_rms_KJL); 
            sprintf('DGO mean = %.3g', mean(pre_err_DGO)); sprintf('DGO RMS = %.3g', pre_rms_DGO); 
            sprintf('ACB mean = %.3g', mean(pre_err_ACB)); sprintf('ACB RMS = %.3g', pre_rms_ACB)}; 
        text(0, 0.5, stext); axis off
    subplot(m, n, ivec(1:end-1) ) 
        scatter(t_KJL, pre_err_KJL, 8 ); hold on; grid on; 
        scatter(t_DGO, pre_err_DGO, 8, 'x'); 
        scatter(t_ACB, pre_err_ACB, 8, 'd'); 
        plot(t_ALL, sigma3, 'g--'); 
        plot(t_ALL, -sigma3, 'g--'); 
        if ~isnan(sigma3) 
            ylim( [ -yrange, yrange ] ); 
        end
        title(ftitle) 
        ylabel(yl)  
        legend('KJL', 'DGO', 'ACB') 

end

function [d_err_STA, d_rms_STA, v_err_STA, v_rms_STA] = calc_res_all(Yobs_STA, Ycalc_STA, DATA)

    % Calculate residuals 
    d_err_STA = Yobs_STA(:,3) - Ycalc_STA(:,3); 
    d_rms_STA = rms(d_err_STA); 
    v_err_STA = Yobs_STA(:,4) - Ycalc_STA(:,4); 
    v_rms_STA = rms(v_err_STA); 

end 


