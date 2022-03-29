

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

% % Atmospheric drag 
% r   = norm(rv0_sat(1:3));            % km 
% const.H   = 88667.0 / 1000;           % m --> km 
% const.r0_drag  = (1000 + const.RN);          % m --> km 
% const.p0  = 3.614e-13 * 1e9;          % kg/m3 --> kg/km^3 
% const.p   = p0*exp( -(r-r0_drag)/H ); 
% const.A   = 15 / 1e6;                 % m^2 --> km^2 
% const.wN  = 2*pi / (16.1*60*60);      % Neptune rotation angular velocity (period = 16.1 hrs)
% 
% const.CD  = 1.88;  

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