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
global muE RE muS AU muM eE wE J2 J3 J4 Cd Cs eop_data 
global A m_SC p p0 r0_drag H CD


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
r   = norm(rv0(1:3));            % km 
const.H   = 88667.0 / 1000;           % m --> km 
const.r0_drag  = (700 + RE);          % m --> km 
const.p0  = 3.614e-13 * 1e9;          % kg/m3 --> kg/km^3 
const.p   = p0*exp( -(r-r0_drag)/H ); 
const.A   = 15 / 1e6;                 % m^2 --> km^2 

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
const.CD  = 1.88;    
   
% ballpark right answer 
rv0  = [   const.RN + 300.25613502827
           const.RN + 0.30084111326
           const.RN + 0.7187226552466
           0.66208631750821
           10.26104887569659
           10.270612922289392 ]; 
   
% ballpark right answer 
rv1  = [   26978.25613502827
          1616.30084111326
          219.7187226552466
         -1.66208631750821
          7.26104887569659
         0.270612922289392 ]; 
   
% ballpark right answer 
rv2  = [   26978.25613502827
          1616.30084111326
          119.7187226552466
         -10.66208631750821
          7.26104887569659
         0.270612922289392 ]; 

nX  = length(rv0); 

%% orbital elements for satellite and Triton 


abcorr  = 'NONE';

%  Convert the epoch to ephemeris time (secs) 
et_t0   = cspice_str2et( t0 );

% get states --> Sun to Earth 
target   = 'Triton';
frame    = 'J2000';
observer = 'Sun';
abcorr   = 'NONE';

% get sun position 
et = et_t0;    % propagate ephemeris time by 1 day in secs 
X_sunE  = spice_state(et, target, frame, abcorr, observer); 

OE0 = rvOrb.rv2orb(rv0, const.muN); 
T0 = 2*pi*sqrt( OE0(1)^3 / const.muN ); 

OE1 = rvOrb.rv2orb(rv1, const.muN); 
T1 = 2*pi*sqrt( OE1(1)^3 / const.muN ); 

OE2 = rvOrb.rv2orb(rv2, const.muN); 
T2 = 2*pi*sqrt( OE2(1)^3 / const.muN ); 


%% integrate EOM 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% Set run state 
disp('Running sim ...')

dt = 10; 

tic
[t, X_NS] = ode45(@fn.EOM, [0 : dt : T0], rv0, options); 
toc 
disp('Pos and Vel end: ')
disp(X_NS(end, 1:6)'); 

% integrate Triton 
abcorr  = 'NONE';

%  Convert the epoch to ephemeris time (secs) 
et_t0   = cspice_str2et( t0 );

% get states --> Neptune to Triton
target   = 'Triton';
frame    = 'J2000';
observer = 'Neptune';
abcorr   = 'NONE';

% get Triton position 
et = et_t0;    % propagate ephemeris time by 1 day in secs 
X_NT0  = spice_state(et, target, frame, abcorr, observer); 

% obtain Triton OE 
OE0_T = rvOrb.rv2orb(X_NT0', const.muT); 

% get angle between satellite and Triton velocities 
r_S = X_NS(1:3); 
r_T = X_NT0(1:3); 
v_S = X_NS(4:6); 
v_T = X_NT0(4:6); 

% angle and velocity angles 
phi_r = acosd( dot(r_S, r_T) / (norm(r_S)*norm(r_T)) ); 
phi_v = acosd( dot(v_S, v_T) / (norm(v_S)*norm(v_T)) ); 

X_NT = []; 
for i = 0 : dt : T0
    
    % propagate by 0.1 day 
    et = et_t0 + i; 
    
    % get velocity 
    X  = spice_state(et, target, frame, abcorr, observer); 
    r_T = X(1:3); 
    v_T = X(4:6); 
    
    X_NT = [X_NT; X]; 
    
end 

% 
% %%
% clear OE 
% for i = 1:length(X_NS)
%     OE(i,:) = fn.rv2orb(X_NS(i,:), const.muN); 
% end 

%% Plot satellite position 

phi_des = 100; 

% [ell_1_min, ell_2_min, amin_AU, emin] = lambert_prob(X_NS, t0, phi_des, plot_option)
[ell_1_min, ell_2_min, amin_AU, emin] = lambert_prob(rv0_T, t0, phi_des, 1)

%% 

% Neptune sphere 
[X, Y, Z] = sphere; 
XN = X * const.RN; YN = Y * const.RN; ZN = Z * const.RN; 

% Triton orbit = OE0_T; 


plot_orbit = 1; 
if plot_orbit == 1
    ftitle = 'Orbit around Neptune'; 
    figure('name', ftitle); 
        surf(XN, YN, ZN); alpha 0.2; shading interp; 
        hold on; grid on; 
        plot3(X_NS(:,1), X_NS(:,2), X_NS(:,3)); 
        plot3(X_NT(:,1), X_NT(:,2), X_NT(:,3));
        scatter3(0, 0, 0, 'filled')
        legend('rv0', 'rv0_t'); 
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

