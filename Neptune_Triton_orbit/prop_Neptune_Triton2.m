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

%% Parameters 

global const 
global eop_data 
global p0 r0_drag H


% ------------------------------------------------------------------------
% Constants 

% const.muE = 398600.4415;          % Earth Gravitational Parameter (km^3/s^2) 
% const.muN = 6.836529e15;          % Neptune GM (m^3 / s^2)
const.muN = 6.836529e15 / ( 1000^3 ); % Neptune GM (km^3 / s^2)

% const.RE  = 6378.1363;            % Earth Radius (km)
const.RN = 24622;                   % Neptune radius (km) 

const.muS = 132712440018;           % Sun’s Gravitational Parameter (km^3/s^2)
const.AU  = 149597870.7;            % 1 Astronomical Unit (km)

% const.muM = 4902.800066;          % Moon’s Gravitational Parameter (km^3/s^2)
G = 6.674e-11;                      % gravitational constant (m^3/kg/s^2)
const.G = G / (1000^3);             % gravitational constant (km^3/kg/s^2)
const.mT = 2.1390e20;               % Triton mass (kg) 
const.muT = const.G * const.mT;     % Triton gravitational parameter (km^3/s^2)  

% const.eE  = 0.081819221456;         % Earth’s eccentricity (orbit? sphere?) 
const.eN = 0.009;                   % Neptune eccentricity (orbit? sphere?) 

% const.wE  = 7.292115146706979e-5;   % Earth’s rotational velocity (rad/s)
% Neptune sidereal day = 16 h, 6 min, 36 s 
const.dayN = (16*60*60) + (6*60) + 36; 
const.wN = 2*pi / const.dayN;       % Neptune angular velocity (rad/s) 

const.m_SC = 2000;                  % satellite mass (kg) 
const.Cd  = 0.04;                   % diffuse reflection 
const.Cs  = 0.04;                   % specular reflection 

eop_data = load('finals_iau1980.txt'); 

% Atmospheric drag 
r   = norm(rv0_sat(1:3));            % km 
const.H   = 88667.0 / 1000;           % m --> km 
const.r0_drag  = (1000 + const.RN);          % m --> km 
const.p0  = 3.614e-13 * 1e9;          % kg/m3 --> kg/km^3 
const.p   = p0*exp( -(r-r0_drag)/H ); 
const.A   = 15 / 1e6;                 % m^2 --> km^2 
const.wN  = 2*pi / (16.1*60*60);      % Neptune rotation angular velocity (period = 16.1 hrs)

const.CD  = 1.88;  

global Cnm Snm 

% Gravity 
Cnm = zeros(181,181);
Snm = zeros(181,181);
fid = fopen('GGM03S.txt','r');
for n=0:180
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);        
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end

%% initial state 

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

% get Triton position 
et = et_t0;    
X_NT0  = spice_state(et, target, frame, abcorr, observer); 
OE0_T = rvOrb.rv2orb(X_NT0', const.muN); 

% obtain desired satellite OE 
a_sat = 1/2 * ( OE0_T(1) + const.RN + 150 ); 
e_sat = 0.9; 
i_sat = OE0_T(3); 
% i_sat = 1.92496219921819; 
argp_sat = OE0_T(4); 
RAAN_sat = OE0_T(5); 
nu_sat = 0; 

% SATELLITE orbital elements and T0 
OE0_sat = [   
        a_sat 
        e_sat 
        i_sat 
        1.16060115269676 
        3.63912435712042 
        0.697239974417234]; 
% OE0_sat = [ 
%         a_sat 
%         e_sat 
%         i_sat 
%         argp_sat 
%         RAAN_sat 
%         nu_sat 
%         ]; 

OE0_test = OE0_T; 
OE0_test(1) = 1/2 * ( OE0_T(1) + const.RN + 150 ); 
OE0_test(2) = 0.9; 
% OE0_test(5) = OE0_sat(5); 

% hyperbolic orbit 
% rv0_test = [19813.3, -16908.2, 2612.7, -22.953, -13.324, 11.316]; 

% "beginning" of hyperbolic orbit 
rv0_test = [          
          1977638.89913003
          2547702.86987409
         -1584975.29348782
         -9.58620042665604
         -12.6791987521143
          7.82653873130335 ]; 
      
      

rv0_sat = rvOrb.orb2rv(OE0_sat, const.muN); 
% rv0_test = rvOrb.orb2rv(OE0_test, const.muN); 
nX  = length(rv0_sat); 

% satellite period 
T_sat = 2*pi*sqrt( OE0_sat(1)^3 / const.muN ); 
T_test = 2*pi*sqrt( OE0_test(1)^3 / const.muN ); 


%% integrate EOM - satellite 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% Set run state 
disp('Running sim ...')

dt = 10; 

tic
[t, X_Nsat] = ode45(@fn.EOM, [0 : dt : T_sat], rv0_sat, options); 
toc 
disp('Pos and Vel end: ')
disp(X_Nsat(end, 1:6)'); 
[t, X_Ntest] = ode45(@fn.EOM, [-T_test : dt : T_test], rv0_test, options); 

% compute rnorm 
for i = 1:length(X_Ntest)
    rnorm_Ntest(i,:) = norm(X_Ntest(i,1:3)); 
end 

% find index of min norm 
i_min = find(rnorm_Ntest == min(rnorm_Ntest)); 
X_Ntest_min = X_Ntest(i_min, :); 
OE_Ntest_min = rvOrb.rv2orb(X_Ntest_min, const.muN); 
OE_Ntest_min(1) = OE0_T(1); 
OE_Ntest_min(2) = 0.9; 
T_test2 = 2*pi*sqrt( OE_Ntest_min(1)^3 / const.muN ); 

rv0_test2 = rvOrb.orb2rv(OE_Ntest_min, const.muN); 
[t, X_Ntest2] = ode45(@fn.EOM, [0: dt : T_test2], rv0_test2, options); 

%% propagate Triton 

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

%%

    labels = {'a', 'e', 'i', 'w (arg of perigee)', 'O (RAAN)', 'nu (true anomaly)'}; 
    ftitle = 'Triton Orbital Elements'; 
    figure('name', ftitle) 
        for i = 1:6 
            subplot(6,1,i) 
            plot([0:dt:T0_T], OE_T(:,i))
            title(labels{i})
        end 
        xlabel('time (s)') 

%% Plot satellite position 

phi_des = 100; 

% [ell_1_min, ell_2_min, amin_AU, emin] = lambert_prob(X_NS, t0, phi_des, plot_option)
% [ell_1_min, ell_2_min, amin_AU, emin] = lambert_prob(rv0_T, t0, phi_des, 1)

%% 

% Neptune sphere 
[X, Y, Z] = sphere; 
XN = X * const.RN; YN = Y * const.RN; ZN = Z * const.RN; 


% Triton orbit = OE0_T; 
rnorm_T = norm(X_NT0(1:3)); 
for theta = 0 : 1 : 360
    X_ecl(theta+1, 1) = cosd(theta) * rnorm_T; 
    X_ecl(theta+1, 2) = sind(theta) * rnorm_T; 
    X_ecl(theta+1, 3) = 0; 
end 

% vernal equinox 
v_eq = [rnorm_T, 0, 0]; 


plot_orbit = 1; 
if plot_orbit == 1
    ftitle = 'Orbit around Neptune'; 
    figure('name', ftitle); 
%         surf(XN, YN, ZN); alpha 0.2; shading interp; 
        ellipsoid(0, 0, 0, const.RN + 1000, const.RN + 1000, const.RN + 1000); alpha 0.2; shading interp; 
        hold on; grid on; axis equal 
%         plot3(X_Nsat(:,1), X_Nsat(:,2), X_Nsat(:,3)); 
        plot3(X_NT(:,1), X_NT(:,2), X_NT(:,3));
%         plot3(X_ecl(:,1), X_ecl(:,2), X_ecl(:,3)); 
%         patch([1 -1 -1 1]*rnorm_T, [1 1 -1 -1]*rnorm_T, [0 0 0 0]*rnorm_T, [1 1 -1 -1]*rnorm_T); alpha 0.2  

        % J2000 axes 
        quiver3(0, 0, 0, v_eq(1), v_eq(2), v_eq(3))
        quiver3(0, 0, 0, 0, rnorm_T, 0); 
        quiver3(0, 0, 0, 0, 0, rnorm_T); 
            txt = 'J2000_x';
            text(v_eq(1),v_eq(2),v_eq(3), txt)
        
        ellipsoid(0, 0, 0, const.RN, const.RN, const.RN); alpha 0.4; shading interp; 
        plot3(X_Ntest(1,1), X_Ntest(1,2), X_Ntest(1,3), 'mo'); 
        plot3(X_Ntest(:,1), X_Ntest(:,2), X_Ntest(:,3), 'm', 'linewidth', 2); 
        plot3(X_Ntest(end,1), X_Ntest(end,2), X_Ntest(end,3), 'm^'); 
        plot3(X_Ntest_min(1), X_Ntest_min(2), X_Ntest_min(3), 'mp'); 

%         plot3(19813.3, -16908.2, 2612.7, 'p')
        plot3(X_Ntest2(1,1), X_Ntest2(1,2), X_Ntest2(1,3), 'co'); 
        plot3(X_Ntest2(:,1), X_Ntest2(:,2), X_Ntest2(:,3), 'c', 'linewidth', 2); 
        plot3(X_Ntest2(end,1), X_Ntest2(end,2), X_Ntest2(end,3), 'c^'); 
%         legend('Neptune', 'sat', 'Triton', 'ecliptic', 'J2000_X', 'J2000_Y', 'J2000_Z'); 
        xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)'); 
        title(ftitle)
        
%     labels = {'a', 'e', 'i', 'w (arg of perigee)', 'O (RAAN)', 'nu (true anomaly)'}; 
%     ftitle = 'Orbital Elements'; 
%     figure('name', ftitle) 
%         for i = 1:6 
%             subplot(6,1,i) 
%             plot(t, OE(:,i))
%             title(labels{i})
%         end 
%         xlabel('time (s)') 
    
end

