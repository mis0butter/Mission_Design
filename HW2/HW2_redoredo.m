%% HW 2 
% Junette Hsin 

% Keplerian elements 
% a     = semi-major axis 
% e     = eccentricity 
% i     = inclination 
% L     = mean longitude 
% wbar  = longitude of perihelion 
% Omega = longitude of ascending node 

close all; clear; 

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

%% transfer angle = 75 deg 

addpath(genpath('mice')); 
addpath(genpath('spice_data')); 

%  Load kernel file 
cspice_furnsh( 'spice_data/naif0011.tls' )
cspice_furnsh( 'spice_data/de421.bsp' )       
cspice_furnsh( 'spice_data/pck00010.tpc ') 

%  Define parameters for a state lookup:
% t0      = 'Oct 20, 2020 11:00 AM CST'; 
t0      = 'May 22, 2000'; 
abcorr  = 'NONE';

%  Convert the epoch to ephemeris time (secs) 
et_t0   = cspice_str2et( t0 );

% get states --> Sun to Earth 
target   = 'Earth';
frame    = 'J2000';
observer = 'Sun';
abcorr   = 'NONE';

% get sun position 
et = et_t0;    % propagate ephemeris time by 1 day in secs 
X_sunE  = spice_state(et, target, frame, abcorr, observer); 

% get states --> Sun to Mars
target   = 'Mars';
frame    = 'J2000';
observer = 'Sun';
abcorr   = 'NONE';

% get sun position 
et = et_t0;    % propagate ephemeris time by 1 day in secs 
X_sunM  = spice_state(et, target, frame, abcorr, observer); 

% get angle between Earth and Mars velocities 
r_E = X_sunE(1:3); 
r_M = X_sunM(1:3); 
v_E = X_sunE(4:6); 
v_M = X_sunM(4:6); 

% angle and velocity angles 
phi_r = acosd( dot(r_E, r_M) / (norm(r_E)*norm(r_M)) ); 
phi_v = acosd( dot(v_E, v_M) / (norm(v_E)*norm(v_M)) ); 

i = 0; 
while abs(phi_r - 75) > 0.01
    
    % propagate by 0.1 day 
    i  = i + 0.01; 
    et = et_t0 + i*86400; 
    
    % get velocity 
    X_sunM  = spice_state(et, target, frame, abcorr, observer); 
    r_M = X_sunM(1:3); 
    v_M = X_sunM(4:6); 

    % get angle 
    phi_r = acosd( dot(r_E, r_M) / (norm(r_E)*norm(r_M)) ); 
    phi_v = acosd( dot(v_E, v_M) / (norm(v_E)*norm(v_M)) ); 
    
end 

% units in AU 
km2AU = 149598073; 
rd = r_E / km2AU; 
ra = r_M / km2AU; 
vd_E = v_E / km2AU; 
va_M = v_M / km2AU; 

    
%% Part 1a
% Lambert - from VALLADO ed. 4, pgs 467 - 475 

% Arrival (Mars), AU units 
ra_mag = norm(ra); 
rd_mag = norm(rd); 

% Vallado method ... 
cos_dv = dot(rd, ra) / (rd_mag * ra_mag); 

% chord 
c = sqrt( rd_mag^2 + ra_mag^2 - 2*rd_mag*ra_mag*cos_dv ); 

% semiperimeter 
s = ( ra_mag + rd_mag + c ) / 2; 

% min semimajor axis 
amin = s/2;  

% initialize 
a_hist  = []; 
dt_hist = []; 
vd_s_hist = []; 
va_s_hist = []; 
vd_l_hist = []; 
va_l_hist = []; 
phi_ds_hist = []; 
phi_dl_hist = []; 
phi_as_hist = []; 
phi_al_hist = []; 

% loop 
for a = amin : 0.01 : 2 

    % time of flight 
    [dt, vd_s, va_s, vd_l, va_l] = a2tof(s, c, a, rd, ra); 
    
    % departure angle, short and long  
    phi_ds = acosd( dot(vd_s, vd_E) / ( norm(vd_s)*norm(vd_E) ) ); 
    phi_dl = acosd( dot(vd_l, vd_E) / ( norm(vd_l)*norm(vd_E) ) ); 

    % arrival angle, short and long 
    phi_as = acosd( dot(va_s, va_M) / ( norm(va_s)*norm(va_M) ) ); 
    phi_al = acosd( dot(va_l, va_M) / ( norm(va_l)*norm(va_M) ) ); 
    
    % save 
    a_hist  = [a_hist; a]; 
    dt_hist = [dt_hist; dt]; 
    vd_s_hist = [vd_s_hist; vd_s]; 
    va_s_hist = [va_s_hist; va_s]; 
    vd_l_hist = [vd_l_hist; vd_l]; 
    va_l_hist = [va_l_hist; va_l]; 
    phi_ds_hist = [phi_ds_hist; phi_ds]; 
    phi_dl_hist = [phi_dl_hist; phi_dl]; 
    phi_as_hist = [phi_as_hist; phi_as]; 
    phi_al_hist = [phi_al_hist; phi_al]; 
    
end 

dt_hist_years = dt_hist / 365; 

% plot 
figure()
    subplot(2,1,1)
        plot(a_hist, dt_hist_years); 
        legend('dt1 short', 'dt1 long', 'dt2 short', 'dt2 long', 'location', 'best'); 
        xlabel('a (AU)'); 
        ylabel('TOF (years)') 
        title('a vs. TOF') 
    subplot(2,1,2) 
        plot(a_hist, phi_ds_hist); hold on; 
        plot(a_hist, phi_dl_hist); 
        plot(a_hist, phi_as_hist); 
        plot(a_hist, phi_al_hist); 
        legend('dep short', 'dep long', 'arr short', 'arr long', 'location', 'best')
        title('a vs. departure/arrival angles'); 
        ylabel('deg') 

% [V1, V2] = LAMBERTBATTIN(rd, ra, 'retro', tof); 

%% Part 1b 
% minimum energy transfer orbit 

pmin = rd_mag*ra_mag/c * (1 - cosd(phi_r)); 

emin = sqrt( 1 - 2*pmin/s ); 

sprintf('a_min = %.5g AU, e_min = %.5g', amin, emin)


%% subfunctions 

function [dt, vd_s, va_s, vd_l, va_l] = a2tof(s, c, a, rd, ra)
% ------------------------------------------------------------------------
% Inputs: 
%   s = semiperimeter (AU) 
%   c = chord (AU) 
%   a = semimajor axis (AU) 
% 
% Outputs (in dt, units days): 
%   dt1_s = time of flight for alpha/beta 1, short
%   dt1_l = time of flight for alpha/beta 1, long 
%   dt2_s = time of flight for alpha/beta 2, short
%   dt2_l = time of flight for alpha/beta 2, long 
% ------------------------------------------------------------------------

% time of flight 
alpha1 = 2 * asind( sqrt( s/(2*a) ) ); 
alpha2 = 360 - alpha1; 
beta1  = 2 * asind( sqrt( (s-c)/(2*a) ) ); 
beta2  = - beta1; 

% alpha1 = acosd( 1 - (norm(rd)+norm(ra)+c)/(2*a) ); 
% alpha2 = 360 - alpha1; 
% beta1  = acosd( 1 - (norm(rd)+norm(ra)-c)/(2*a) ); 
% beta2  = - beta1; 

% sun mu (m^3/s^2)
mu_sun_m3s2 = 1.32712440018e20; 
m2au = 6.68459e-12; % au per m 
mu_sun_aus = mu_sun_m3s2 * m2au^3; 

% time of flight 
nrev = 0; 
dt1_s = sqrt( a^3/mu_sun_aus ) * ... 
    (( 2*nrev*180 + alpha1 - sind(alpha1) - (beta1 - sind(beta1)) )); 
% dt1_l = sqrt( a^3/mu_sun_aus ) * ... 
%     (( 2*nrev*180 + alpha1 - sind(alpha1) + (beta1 - sind(beta1)) )); 
% dt2_s = sqrt( a^3/mu_sun_aus ) * ... 
%     (( 2*nrev*180 + alpha2 - sind(alpha2) - (beta2 - sind(beta2)) )); 
dt2_l = sqrt( a^3/mu_sun_aus ) * ... 
    (( 2*nrev*180 + alpha2 - sind(alpha2) + (beta2 - sind(beta2)) )); 

% dt = [dt1_s, dt1_l, dt2_s, dt2_l]; 
dt = [dt1_s, dt2_l]; 

% convert from seconds to days 
dt = dt / (60*60*24); 

% velocities 
DE1   = alpha1 - beta1;
F1    = 1.0 - (a/norm(rd)) * (1.0 - cosd(DE1) );
GDot1 = 1.0 - (a/norm(ra)) * (1.0 - cosd(DE1) );
G1_s  = dt1_s - sqrt(a*a*a/mu_sun_aus) * (DE1 - sind(DE1));
for i = 1:3
    vd1_s(i) = ( ra(i) - F1*rd(i) )/G1_s;
    va1_s(i) = ( GDot1*ra(i) - rd(i) )/G1_s;
end
% G1_l  = dt1_l - sqrt(a*a*a/mu_sun_aus) * (DE1 - sind(DE1));
% for i = 1:3
%     vd1_l(i) = ( ra(i) - F1*rd(i) )/G1_l;
%     va1_l(i) = ( GDot1*ra(i) - rd(i) )/G1_l;
% end

DE2   = alpha2 - beta2; 
F2    = 1.0 - (a/norm(rd)) * (1.0 - cosd(DE2) ); 
GDot2 = 1.0 - (a/norm(ra)) * (1.0 - cosd(DE2) ); 
% G2_s  = dt2_s - sqrt(a*a*a/mu_sun_aus) * (DE2 - sind(DE2)); 
% for i = 1:3
%     vd2_s(i) = ( ra(i) - F2*rd(i) )/G2_s; 
%     va2_s(i) = ( GDot2*ra(i) - rd(i) )/G2_s; 
% end 
G2_l  = dt2_l - sqrt(a*a*a/mu_sun_aus) * (DE2 - sind(DE2)); 
for i = 1:3
    vd2_l(i) = ( ra(i) - F2*rd(i) )/G2_l; 
    va2_l(i) = ( GDot2*ra(i) - rd(i) )/G2_l; 
end 

%% VELOCITIES 

rd_mag = norm(rd); 
ra_mag = norm(ra); 
cos_dv = dot(ra, rd) / ( rd_mag*ra_mag );
dv     = acosd(cos_dv); 
sin_dv = dt1_s * sqrt( 1 - cos_dv^2 ); 

l = (rd_mag+ra_mag) / ( 4*sqrt( rd_mag*ra_mag ) * cosd(dv/2) ) - 1/2; 
m = mu_sun_aus * dt1_s^2 / ( 2*sqrt( rd_mag*ra_mag ) * cosd(dv/2) )^3;

% LOOP
y  = 1; 
x1 = m/y^2 - l; 
x2 = 4/3 * ( 1 + 6*x1/5 + 6*8*x1^2/(5*7) + 6*8*10*x1^3/(5*7*9) ); 

[vd_s, va_s] = LAMBERTBATTIN_au(rd, ra, 'pro', dt1_s); 
[vd_l, va_l] = LAMBERTBATTIN_au(rd, ra, 'pro', dt2_l); 



end 


function [r_ecl] = xyz_ecl(Teph, planet)
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

end 


%% Kepler equation solver 

function E = keplerEq(M,e,eps)
% Function solves Kepler's equation M = E-e*sin(E)
% Input - Mean anomaly M [rad] , Eccentricity e and Epsilon 
% Output  eccentric anomaly E [DEG]. 
   	En  = M;
	Ens = En - (En-e*sind(En)- M)/(1 - e*cosd(En));
	while ( abs(Ens-En) > eps )
		En = Ens;
		Ens = En - (En - e*sind(En) - M)/(1 - e*cosd(En));
    end
	E = Ens;
end

function E = kepler(M, e)
    f = @(E) E - e * sind(E) - M;
    E = fzero(f, M);  % <-- I would use M as the initial guess instead of 0
end

% 

function h2_text(h2, text1, text2, text3, text4) 

    pos = get(h2, 'position');     
    delete(findall(gcf,'type','annotation')); 
    text = { ''; ''; text1; text2; text3; text4 }; 
    annotation('textbox', pos, ...
      'String', text, ...
      'edgecolor', 'none');
    axis off 

end 


function rv = spice_state(epoch, target, frame, abcorr, observer) 

    rv = zeros(length(epoch), 6); 
    
    for i = 1:length(epoch) 

        %  Look-up the state for the defined parameters.
        starg   = mice_spkezr( target, epoch(i), frame, abcorr, observer);
        rv(i,:) = starg.state(1:6); 
        
    end 

end 


