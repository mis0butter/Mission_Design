% HW 5 
% Junette Hsin 

% finalproj_EOM_batch; 
close all; format long g; 
load WS_batch.mat 

% Load SPICE stuff 
addpath(genpath('mice')); 
addpath(genpath('spice_data')); 
cspice_furnsh( 'spice_data/naif0011.tls' )
cspice_furnsh( 'spice_data/de421.bsp' )       
cspice_furnsh( 'spice_data/pck00010.tpc ' ) 

%% 

global eopdata const

SAT_Const

% read Earth orientation parameters
fid = fopen('eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

%%

% plot options 
plot_orbit = 0; 
plot_RSW   = 0; 

% Selct sim params 
run_STA  = 1 ;      % run STATIONS 
run_DATA = 0 ;      % run DATA 
save_nag = 0 ;      % save NAG 
DAYS     = 1;       % PROPAGATE to 1 day, 3 days, or 7 days. 

%% EKF - all observations 

% loop through STATIONS 
if run_STA == 1
    DATA = 0;       % 0 = all data 
    for STATIONS = 3    % 0 = all stations, 1 = KJL, 2 = DGO, 3 = ACB
        EKF_STA_DATA; 
    end
end 

% loop through DATA 
if run_DATA == 1
    STATIONS = 0;       % 0 = all stations 
    for DATA = 1:2      % 0 = all data, 1 = range, 2 = range-rate, 3 = short-arc 
        EKF_STA_DATA; 
    end
end

X_sol1 = [ -6330.16736001325
            3306.81591178162
            127.736863438565
           -3.43787907681733
           -6.63350511630163
           -0.235613730204275 ]; 


%% save NAG 

% Save in NAG format 
if save_nag == 1
    
    load WS_R_7.mat     % case A 
        save_case(Xi, P_bar, 'A'); 
    load WS_RR_7.mat    % case B 
        save_case(Xi, P_bar, 'B'); 
    load WS_KJL_7.mat   % case C 
        save_case(Xi, P_bar, 'C'); 
    load WS_DGO_7.mat   % case D 
        save_case(Xi, P_bar, 'D'); 
    load WS_ACB_7.mat   % case E 
        save_case(Xi, P_bar, 'E'); 
    load WS_ALL_7.mat   % case F (LONG arc now) 
        save_case(Xi, P_bar, 'F'); 
    load WS_ALL_7_short_arc.mat   % case F (short arc now) 
        save_case(Xi, P_bar, 'G'); 
    
    save('hsinj.mat' ... 
        , 'hsinj_pos_caseA', 'hsinj_poscov_caseA' ...   % Range only 
        , 'hsinj_pos_caseB', 'hsinj_poscov_caseB' ...   % Range-rate only 
        , 'hsinj_pos_caseC', 'hsinj_poscov_caseC' ...   % Kwajalein only 
        , 'hsinj_pos_caseD', 'hsinj_poscov_caseD' ...   % Diego-Garcia only 
        , 'hsinj_pos_caseE', 'hsinj_poscov_caseE' ...   % Arecibo only 
        , 'hsinj_pos_caseF', 'hsinj_poscov_caseF' ...   % All stations, all data 
        , 'hsinj_pos_caseG', 'hsinj_poscov_caseG' ... 
    );  
    
end

%% subfunctions 

function save_case(Xi, P_bar, case_str)

pos_str = ['hsinj_pos_case' case_str]; 
poscov_str = ['hsinj_poscov_case' case_str]; 

evalin('base', [pos_str, ' = Xi(1:3)'])     % Not recommended but oh well 
evalin('base', [poscov_str, ' = P_bar(1:3,1:3)'])     % Not recommended but oh well 

end








