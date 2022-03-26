clc
clear
format long g

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

year = 1980;
month = 10;
day = 2;
hour = 23;
minute = 41;
second = 24.1137599987633;
MJD_UTC = Mjday(year, month, day, hour, minute, second);
rv_eci = 1e3*[2769.30588515169, -6067.79991631, -321.35732976219, 1.95719428899467, 1.21970228570168, -7.35245769612438]; % [m]
rv_ecef = ECI2ECEF(MJD_UTC, rv_eci)
rv_eci = ECEF2ECI(MJD_UTC, rv_ecef')

%% HW 4

rv_ECEF = [ -28738.3218400000; -30844.0723200000; -6.71800000000000; 0; 0; 0 ];
r_ECEF = [ -28738.3218400000; -30844.0723200000; -6.71800000000000 ]; 
JD = 2458088.50055556; 
MJD = JD - 2400000.5; 

%%

eop_data = load('finals_iau1980.txt'); 

[r_ECI] = fn.ECEFtoECI(eop_data, JD, r_ECEF); 
%%

X_KJL_ECI  = []; 
X_KJL_ECEF = [-6143584  1364250  1033743 0 0 0] / 1000;  % Kwajalein 

i = 0; 
for t = MJD : 60 : MJD + 3600 
    i = i + 1; 
    X_KJL_ECI(:,i) = ECEF2ECI(t, X_KJL_ECEF); 
end





