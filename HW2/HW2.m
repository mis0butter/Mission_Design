%% HW 2 
% Junette Hsin 

% Keplerian elements 
% a     = semi-major axis 
% e     = eccentricity 
% i     = inclination 
% L     = mean longitude 
% wbar  = longitude of perihelion 
% Omega = longitude of ascending node 

clear; clc 

% Earth 
Earth.a0 =      1.00000261; 
Earth.da =      0.00000562; 
Earth.e0 =      0.01671123; 
Earth.de =     -0.00004392; 
Earth.I0 =     -0.00001531; 
Earth.dI =     -0.01294668;
Earth.L0 =    100.46457166; 
Earth.dL =  35999.37244981; 
Earth.wbar0 = 102.93768193; 
Earth.dwbar =   0.32327364;
Earth.Omega0 =  0; 
Earth.dOmega =  0; 

% Mars 
Mars.a0 =      1.52371034; 
Mars.da =      0.00001847; 
Mars.e0 =      0.09339410; 
Mars.de =      0.00007882; 
Mars.I0 =      1.84969142; 
Mars.dI =     -0.00813131; 
Mars.L0 =     -4.55343205; 
Mars.dL =  19140.30268499; 
Mars.wbar0 = -23.94362959; 
Mars.dwbar =   0.44441088; 
Mars.Omega0 = 49.55953891; 
Mars.dOmega = -0.29257343; 

%% Try Future Earth 

% Departure (Earth), AU units 
Teph = 2451545.0;
rd = xyz_ecl(Teph, Earth); 
rd_mag = norm(rd); 

% transfer angle (deg)
phi = 75; 

% Arrival (Mars), AU units 
tof = 26*7*24*60*60; 
Teph = 2451545.0 + tof;   % 26 weeks 
ra = xyz_ecl(Teph, Mars); 
ra_mag = norm(ra); 

% Vallado method ... 
cos_dv = dot(rd, ra); 
cos_dv = cos_dv / norm(cos_dv); 

% chord 
c = sqrt( rd_mag^2 + ra_mag^2 - 2*rd_mag*ra_mag*cos_dv ); 

% semiperimeter 
s = ( ra_mag + rd_mag + c ) / 2; 

% min semimajor axis 
a = s/2;  

% time of flight 
ae = 2 * asind( sqrt( s/(2*a) ) ); 
be = 2 * asind( sqrt( (s-c)/(2*a) ) ); 
% dt = sqrt( a^3/mu ) * (( 2*nrev*pi + ae - sin(ae) + (be - sin(be)) )); 

% [V1, V2] = LAMBERTBATTIN(rd, ra, 'retro', tof); 




%% Try Future Earth 

function [r_ecl] = xyz_ecl(Teph, planet)

T = (Teph - 2451545.0)/36525;

% propagate 6 elements 
a = planet.a0 + planet.da * T; 
e = planet.e0 + planet.de * T; 
I = planet.I0 + planet.dI * T; 
L = planet.L0 + planet.dL * T; 
wbar = planet.wbar0 + planet.dwbar * T; 
Omega = planet.Omega0 + planet.dOmega * T; 

% argument of perihelion 
w = planet.wbar0 - planet.Omega0; 

% mean anomaly 
M = planet.L0 - planet.wbar0; 

% mean motion 
n = 2*pi / 2; 

% eccentric anomaly 
E = keplerEq(M, planet.e0, eps); 

% heliocentric coordinates in orbital plane, xp algined from focus to
% perihelion 
xp = a * (cos(E) - e); 
yp = a * sqrt( 1 - e^2 ) * sin(E); 
zp = 0; 

% coordinates in J2000 ecliptic plane, x aligned toward equinox 
x_ecl = (  cos(w)*cos(Omega) - sin(w)*sin(Omega)*cos(I) ) * xp + ... 
        ( -sin(w)*cos(Omega) - cos(w)*sin(Omega)*cos(I) ) * yp ; 
y_ecl = (  cos(w)*sin(Omega) + sin(w)*cos(Omega)*cos(I) ) * xp + ... 
        ( -sin(w)*sin(Omega) + cos(w)*cos(Omega)*cos(I) ) * yp;     
z_ecl = ( sin(w)*sin(I) ) * xp + ( cos(w)*sin(I) ) * yp; 
r_ecl = [x_ecl; y_ecl; z_ecl]; 

% obliquity at J2000 (deg) 
e_obl = 23.43928; 

% "ICRF" or "J2000 frame"
x_eq  = x_ecl; 
y_eq  = cos(e_obl)*y_ecl - sin(e_obl)*z_ecl; 
z_eq  = sin(e_obl)*y_ecl + cos(e_obl)*z_ecl; 

end 


%% Kepler equation solver 

function E = keplerEq(M,e,eps)
% Function solves Kepler's equation M = E-e*sin(E)
% Input - Mean anomaly M [rad] , Eccentricity e and Epsilon 
% Output  eccentric anomaly E [rad]. 
   	En  = M;
	Ens = En - (En-e*sin(En)- M)/(1 - e*cos(En));
	while ( abs(Ens-En) > eps )
		En = Ens;
		Ens = En - (En - e*sin(En) - M)/(1 - e*cos(En));
    end
	E = Ens;
end






