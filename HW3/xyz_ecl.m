function [r_ecl, r_eq, oe] = xyz_ecl(Teph, planet)
% ------------------------------------------------------------------------
% Teph  = ephemeris time (days) 
% r_ecl = J2000 ecliptic plane coordinates 
% ------------------------------------------------------------------------

T = (Teph - 2451545.0)/36525;

% propagate 6 elements 
a = planet.a0 + planet.da * T; 
e = planet.e0 + planet.de * T; 
I = planet.I0 + planet.dI * T; 
L = planet.L0 + planet.dL * T; 
wbar = planet.wbar0 + planet.dwbar * T; 
Omega = planet.Omega0 + planet.dOmega * T; 

oe = [a; e; I; L; wbar; Omega]; 

% argument of perihelion 
w = wbar - Omega; 

% mean anomaly 
M = L - wbar; 
if abs(M) > 180
    if M > 0 
        rem = mod(M, 180); 
        M = - 180 + rem; 
    else
        rem = mod(-M, 180); 
        M = 180 - rem; 
    end 
end 

% eccentric anomaly 
% E = keplerEq(M, e, eps); 
E = kepler(M, 180/pi*e); 

% heliocentric coordinates in orbital plane, xp algined from focus to
% perihelion 
xp = a * (cosd(E) - e); 
yp = a * sqrt( 1 - e^2 ) * sind(E); 
zp = 0; 

% coordinates in J2000 ecliptic plane, x aligned toward equinox 
x_ecl = (  cosd(w)*cosd(Omega) - sind(w)*sind(Omega)*cosd(I) ) * xp + ... 
        ( -sind(w)*cosd(Omega) - cosd(w)*sind(Omega)*cosd(I) ) * yp ; 
y_ecl = (  cosd(w)*sind(Omega) + sind(w)*cosd(Omega)*cosd(I) ) * xp + ... 
        ( -sind(w)*sind(Omega) + cosd(w)*cosd(Omega)*cosd(I) ) * yp;     
z_ecl = ( sind(w)*sind(I) ) * xp + ( cosd(w)*sind(I) ) * yp; 
r_ecl = [ x_ecl; y_ecl; z_ecl ]; 

% obliquity at J2000 (deg) 
e_obl = 23.43928; 

% "ICRF" or "J2000 frame"
x_eq  = x_ecl; 
y_eq  = cosd(e_obl)*y_ecl - sind(e_obl)*z_ecl; 
z_eq  = sind(e_obl)*y_ecl + cosd(e_obl)*z_ecl; 

r_eq = [x_eq; y_eq; z_eq]; 

end 

function E = kepler(M, e)
    f = @(E) E - e * sind(E) - M;
    E = fzero(f, M);  % <-- I would use M as the initial guess instead of 0
end