%% Junette Hsin 

% 0 = no plot. 1 = plot animation. 2 = plot all at once 
plot_option = 2; 
 
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


%% Problem 1 

% (Ground Tracks) Consider the International Space Station (look up its orbital elements)
% in a Keplerian orbit. Assume that at t = 0, the ISS, its orbital ascending node, and the
% Greenwich Meridian are coincident. Note that all references to node crossings and ground
% tracks are from a viewpoint attached to the surface of the Earth.

% Draw and compare the ground tracks for 15-minute duration in two cases: 
% i. A non-rotating Earth
% ii. A uniformly rotating Earth.

% t = 0, ascending node and Greenwich Meridian coincident 

% ISS OEs (from https://in-the-sky.org/spacecraft_elements.php?id=25544)  
e0 = 0.00048;   
i0 = 51.644 * pi/180;    % deg --> rad
w0 = 30.4757 * pi/180;   % deg --> rad
% w0 = 0; 
O0 = 0 * pi/180;         % deg --> rad
% M0 = 39.7178 * pi/180;   % deg --> rad (should be true anomaly) 
M0 = (2*pi - w0);        % rad (should be true anomaly) 

% mean motion --> semimajor axis 
n = 15.50094 / 86400 * (2*pi);   % rev/day --> rad/s 
mu_E_m3 = 3.986004418e14;  % m^3/s^2
mu_E_km3 = mu_E_m3 * ( 1e-3 )^3 ; 
a0 = (mu_E_km3 / n^2)^(1/3); 

R_E = 6378.1370;        % km 
w_E = 7.292115e-5;      % rad/s 

oe0 = [a0; e0; i0; w0; O0; M0]; 
rv0 = rvOrb.orb2rv(oe0, mu_E_km3); 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% propagate orbit 
[t, rv_a] = ode45(@fn.EOM, [0 : 15*60], rv0, options); 


% ------------------------------------------------------------------------
% problem 1.a.i

lla = []; 
lla_1_a = []; 
for i = 1:length(rv_a)
    
%     lla(i,:) = ecef2lla(rv(i,1:3)); 
    lla_1_a(i,:) = ecef2lla_1(rv_a(i,:)', R_E, mu_E_km3); 
    lla_1_a(i, 1:2) = lla_1_a(i, 1:2) * 180/pi; 
    
end 
% fn.plot3_xyz(rv);

% ------------------------------------------------------------------------
% problem 1.a.ii 

% Axis 3 rotation matrix 
[lla_rot_a, rv_rot_a] = lla_rv_rot(t, rv_a, w_E, R_E, mu_E_km3); 


% ------------------------------------------------------------------------
% PLOT 
fname = '1.a Rotating and Non-Rotating Earth'; 
plot_option = 2; 
plot_gt(plot_option, fname, lla_1_a, lla_rot_a)


%% problem 1.b

i0 = (180 - 51.644) * pi/180;    % deg --> rad

% mean motion --> semimajor axis 
n = 15.50094 / 86400 * (2*pi);   % rev/day --> rad/s 
mu_E_m3 = 3.986004418e14;  % m^3/s^2
mu_E_km3 = mu_E_m3 * ( 1e-3 )^3 ; 
a0 = (mu_E_km3 / n^2)^(1/3); 

R_E = 6378.1370;        % km 
w_E = 7.292115e-5;      % rad/s 

oe0 = [a0; e0; i0; w0; O0; M0]; 
rv0 = rvOrb.orb2rv(oe0, mu_E_km3); 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% propagate orbit 
[t, rv_b] = ode45(@fn.EOM, [0 : 15*60], rv0, options); 


% ------------------------------------------------------------------------
% Problem 1.b.i 

lla = []; 
lla_1_b = []; 
for i = 1:length(rv_b)
    
%     lla_1_b(i,:) = ecef2lla(rv(i,1:3)); 
    lla_1_b(i,:) = ecef2lla_1(rv_b(i,:)', R_E, mu_E_km3); 
    lla_1_b(i, 1:2) = lla_1_b(i, 1:2) * 180/pi; 
    
end 
% fn.plot3_xyz(rv); 


% ------------------------------------------------------------------------
% problem 1.b.ii 

% Axis 3 rotation matrix 
[lla_rot_b, rv_rot_b] = lla_rv_rot(t, rv_b, w_E, R_E, mu_E_km3); 


% ------------------------------------------------------------------------
% PLOT 
fname = '1.b Rotating and Non-Rotating Earth'; 
plot_option = 2; 
plot_gt(plot_option, fname, lla_1_b, lla_rot_b)



%% problem 1.c.i 

% PROGRADE 

% Constants
mu = mu_E_km3;
J2 = 1.082e-3;

% O Precession Calcs
Odot = -(3/2)*n*(R_E / a0)^2 * J2 * (1/(1-e0^2)^(1/2)) * cos(i0); % O precession
sprintf('Odot precession: %.3f', Odot)

% ISS OEs (from https://in-the-sky.org/spacecraft_elements.php?id=25544)  
i0 = 51.644 * pi/180;    % deg --> rad
M0 = 0; 

oe0 = [a0; e0; i0; w0; O0; M0]; 
rv0 = rvOrb.orb2rv(oe0, mu_E_km3); 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% propagate orbit 
T = 2*pi*sqrt(a0^3/mu_E_km3); 
[t, rv_c] = ode45(@fn.EOM_J2, [0 : T], rv0, options); 

% final oe 
oef = rvOrb.rv2orb(rv_c(end,:), mu_E_km3); 
Odot = oe0(5) - oef(5); 

lla = []; 
lla_1_c = []; 
for i = 1:length(rv_c)
    
%     lla(i,:) = ecef2lla(rv(i,1:3)); 
    lla_1_c(i,:) = ecef2lla_1(rv_c(i,:)', R_E, mu_E_km3); 
    lla_1_c(i, 1:2) = lla_1_c(i, 1:2) * 180/pi; 
    
end 
% fn.plot3_xyz(rv); 


% ------------------------------------------------------------------------
% Axis 3 rotation matrix 
[lla_rot_c, rv_rot_c] = lla_rv_rot(t, rv_c, w_E, R_E, mu_E_km3); 


% ------------------------------------------------------------------------
% PLOT

fname = '1.c Fixed and Rotating Earth with Secular Precession Prograde'; 
plot_option = 2; 
plot_gt(plot_option, fname, lla_1_c, lla_rot_c)

% ------------------------------------------------------------------------
% RETROGRADE 

% Constants
mu = mu_E_km3;
J2 = 1.082e-3;

% O Precession Calcs
Odot = -(3/2)*n*(R_E / a0)^2 * J2 * (1/(1-e0^2)^(1/2)) * cos(i0); % O precession
sprintf('Odot precession: %.3f', Odot)

% ISS OEs (from https://in-the-sky.org/spacecraft_elements.php?id=25544)  
i0 = (180 - 51.644) * pi/180;    % deg --> rad

oe0 = [a0; e0; i0; w0; O0; M0]; 
rv0 = rvOrb.orb2rv(oe0, mu_E_km3); 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% propagate orbit 
T = 2*pi*sqrt(a0^3/mu_E_km3); 
[t, rv_c] = ode45(@fn.EOM_J2, [0 : T], rv0, options); 

% final oe 
oef = rvOrb.rv2orb(rv_c(end,:), mu_E_km3); 
Odot = oe0(5) - oef(5); 

lla = []; 
lla_1_c = []; 
for i = 1:length(rv_c)
    
%     lla(i,:) = ecef2lla(rv(i,1:3)); 
    lla_1_c(i,:) = ecef2lla_1(rv_c(i,:)', R_E, mu_E_km3); 
    lla_1_c(i, 1:2) = lla_1_c(i, 1:2) * 180/pi; 
    
end 
% fn.plot3_xyz(rv); 


% ------------------------------------------------------------------------
% Axis 3 rotation matrix 
[lla_rot_c, rv_rot_c] = lla_rv_rot(t, rv_c, w_E, R_E, mu_E_km3); 


% ------------------------------------------------------------------------
% PLOT

fname = '1.c Fixed and Rotating Earth with Secular Precession Retrograde '; 
plot_option = 2; 
plot_gt(plot_option, fname, lla_1_c, lla_rot_c)


%% problem 1.c.ii

% Constants
mu = mu_E_km3;
J2 = 1.082e-3;

% O Precession Calcs
Odot = -(3/2)*n*(R_E / a0)^2 * J2 * (1/(1-e0^2)^(1/2)) * cos(i0); % O precession
sprintf('Odot precession: %.3f', Odot)

% ISS OEs (from https://in-the-sky.org/spacecraft_elements.php?id=25544)  
i0 = 51.644 * pi/180;    % deg --> rad

oe0 = [a0; e0; i0; w0; O0; M0]; 
rv0 = rvOrb.orb2rv(oe0, mu_E_km3); 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% propagate orbit 
T = 2*pi*sqrt(a0^3/mu_E_km3); 
[t, rv_c] = ode45(@fn.EOM_J2, [0 : T], rv0, options); 

lla = []; 
lla_1_c = []; 
for i = 1:length(rv_c)
    
%     lla(i,:) = ecef2lla(rv(i,1:3)); 
    lla_1_c(i,:) = ecef2lla_1(rv_c(i,:)', R_E, mu_E_km3); 
    lla_1_c(i, 1:2) = lla_1_c(i, 1:2) * 180/pi; 
    
end 
% fn.plot3_xyz(rv); 


% ------------------------------------------------------------------------
% Axis 3 rotation matrix 
[lla_rot_c, rv_rot_c] = lla_rv_rot(t, rv_c, w_E, R_E, mu_E_km3); 


% ------------------------------------------------------------------------
% PLOT

fname = '1.c Rotating and Non-Rotating Earth'; 
plot_option = 2; 
plot_gt(plot_option, fname, lla_1_c, lla_rot_c)


% ------------------------------------------------------------------------
% Precession matches up? 
oe_a = rvOrb.rv2orb(rv_a(end,:), mu_E_km3); 
oe_c = rvOrb.rv2orb(rv_c(end,:), mu_E_km3); 
% oe_d = rvOrb.rv2orb(rv_c(901,:), mu_E_km3); 

O_a = oe_a(5); 
O_c = oe_c(5); 



%% Problem 1.c.ii another way  

oe_c0 = rvOrb.rv2orb(rv0, mu_E_km3); 
dt = 1; 

Odot = -3/2 * n * J2 * ( R_E / norm(rv0(1:3)) )^2 * 1/( 1-e0^2 )^2 * cos(i0); 
wdot = -3/4 * n * (R_E / norm(rv0(1:3)))^2 * J2 * (1 - 5*cosd(i0)^2) / (1-e0^2)^2;
Mdot = n * (1-3/4*(R_E / norm(rv0(1:3)))^2 * J2 * (1 - 3*cosd(i0)^2) / (1-e0^2)^(3/2));

oe_c = oe_c0; 
rv_c2 = rv0'; 
for i = 1 : T
    
    oe = oe_c0; 
    
    % augment mean anomaly 
    M0 = oe_c0(6); 
    M = M0 + Mdot * (dt*i); 
    M = mod(M, 2*pi); 
    oe(6) = M; 
    
    % augment RAAN 
    O0 = oe_c0(5); 
    O = O0 + Odot * (dt*i); 
    O = mod(O, 2*pi); 
    oe(5) = O; 
    
    % augment perigee 
    w0 = oe_c0(4); 
    w = w0 + wdot * (dt*i); 
    w = mod(w, 2*pi); 
    oe(4) = w; 
    
    oe_c(i+1,:) = oe; 
    rv_c2(i+1,:) = rvOrb.orb2rv(oe, mu_E_km3); 
    
end 

lla = []; 
lla_1_c2 = []; 
for i = 1:length(rv_c2)
    
%     lla(i,:) = ecef2lla(rv(i,1:3)); 
    lla_1_c2(i,:) = ecef2lla_1(rv_c2(i,:)', R_E, mu_E_km3); 
    lla_1_c2(i, 1:2) = lla_1_c2(i, 1:2) * 180/pi; 
    
end 



% ------------------------------------------------------------------------
% Axis 3 rotation matrix 
[lla_rot_c2, rv_rot_c2] = lla_rv_rot(t, rv_c2, w_E, R_E, mu_E_km3); 
    
for i = 1:length(rv_rot_c2)
%     lla(i,:) = ecef2lla(rv(i,1:3)); 
    lla_rot_c2(i,:) = ecef2lla_1(rv_rot_c2(i,:)', R_E, mu_E_km3); 
    lla_rot_c2(i, 1:2) = lla_rot_c2(i, 1:2) * 180/pi; 
    
end 

% fn.plot3_xyz(rv); 
% ------------------------------------------------------------------------
% PLOT

fname = '1.c Fixed and Rotating Earth with Secular Precession'; 
plot_option = 2; 
plot_gt(plot_option, fname, lla_1_c2, lla_rot_c2)


%% Problem 1.c.ii ANOTHER ANOTHER way 

Td = 2*pi / (wdot + Mdot); 
Dn = 2*pi / ( w_E - Odot ); 

u0 = w0 + M0; 

u = u0 + udot * T + 2*pi/udot; 

dlambda = (-w_E - Odot) * Td; 




%% subfunctions 

function [lla_rot, rv_rot] = lla_rv_rot(t, rv_a, w_E, R_E, mu_E_km3)

    theta0 = 0; 
    dt = t(2) - t(1); 

    lla_rot = ecef2lla_1(rv_a(1,:)', R_E, mu_E_km3)'; 
    lla_rot(1, 1:2) = lla_rot(1, 1:2) * 180/pi; 
    rv_rot = rv_a(1,:); 
    
    for i = 2:length(rv_a)

        theta = -w_E * dt * (i-1) + theta0; 
        r_rot = fn.rotate_xyz(rv_a(i,1:3), theta, 3); 
        v_rot = fn.rotate_xyz(rv_a(i,4:6), theta, 3); 
        rv_rot(i,:) = [r_rot; v_rot]; 

        lla_rot(i,:) = ecef2lla_1(rv_rot(i,:)', R_E, mu_E_km3); 
        lla_rot(i, 1:2) = lla_rot(i, 1:2) * 180/pi; 


    end 

end 

% ------------------------------------------------------------------------
function plot_gt(plot_option, fname, lla_1, lla_rot)

%     fname = '1.a Rotating Earth'; 
if plot_option == 1

    figure('name', fname)
    plot(lla_1(:, 2), lla_1(:, 1), 'b'); 
    hold on; grid on; 
    lh = plot(lla_rot(:, 2), lla_rot(:, 1), 'r'); 
    xl = xlim; yl = ylim; 

    cla;  
    xlim(xl); ylim(yl); 
    for i = 1:length(lla_rot)

        cla 

        % fixed Earth 
        plot(lla_1(1,2), lla_1(1,1), 'bo'); 
        lh1 = plot(lla_1(1:i, 2), lla_1(1:i, 1), 'b'); 
        plot(lla_1(i, 2), lla_1(i, 1), 'b^'); 

        % rotating Earth 
        plot(lla_rot(1,2), lla_rot(1,1), 'ro'); 
        lh2 = plot(lla_rot(1:i, 2), lla_rot(1:i, 1), 'r'); 
        lh = plot(lla_rot(i, 2), lla_rot(i, 1), 'r^'); 

        legend([lh1 lh2], 'fixed E', 'rot E') 

        title(sprintf( '%d / %d', i, length(lla_rot) )); 
        pause(0.001) 

    end 
    
    title(fname) 
    xlabel('Longitude (deg)'); 
    ylabel('Latitude (deg)');

elseif plot_option == 2

    figure('name', fname)
    hold on; grid on; 
    
    % fixed Earth 
    plot(lla_1(1,2), lla_1(1,1), 'bo'); 
    lh1 = plot(lla_1(:, 2), lla_1(:, 1), 'b'); 
    plot(lla_1(end, 2), lla_1(end, 1), 'b^'); 

    % rotating Earth 
    plot(lla_rot(1,2), lla_rot(1,1), 'ro'); 
    lh2 = plot(lla_rot(:, 2), lla_rot(:, 1), 'r'); 
    lh = plot(lla_rot(end, 2), lla_rot(end, 1), 'r^'); 

    legend([lh1 lh2], 'fixed E', 'rot E') 

    title(fname) 
    xlabel('Longitude (deg)'); 
    ylabel('Latitude (deg)');
    
end 

end 

% ------------------------------------------------------------------------
% ECEF to geodetic (lat, lon) 
function lla = ecef2lla_1(rv_ecef, R_E, mu)

% extract data 
r_x = rv_ecef(1); r_y = rv_ecef(2); r_z = rv_ecef(3); 
r_norm = norm(rv_ecef); 

% obtain OEs 
oe = rvOrb.rv2orb(rv_ecef, mu); 
a = oe(1); e = oe(2); i = oe(3); 
w = oe(4); O = oe(5); nu = oe(6); 

% solve for longitude 
r_delta = sqrt( r_x^2 + r_y^2 ); 
a_sin = asin(r_y / r_delta); 
a_cos = acos(r_x / r_delta); 

% find quadrant 
if a_sin > 0
    if a_cos > 0
        a = abs(a_cos); 
    else
        a = pi - abs(a_cos); 
    end 
else
    if a_cos > 0
        a = 2*pi - abs(a_cos); 
    else
        a = pi + abs(a_cos); 
    end     
end 

lon = mod(a, 2*pi); 

% iterate for latitude 
delta = asin(r_z / r_norm);

% assume first guess 
lat0 = delta; 
lat  = lat0; 

C = R_E / sqrt( 1 - e^2 * ( sin(lat) )^2 ); 
S = R_E * (1-e^2) / sqrt( 1 - e^2*( sin(lat) )^2 ); 

tan_lat = (r_z + C * e^2 * sin(delta)) / r_delta; 
lat = atan(tan_lat); 
h = r_z / ( sin(delta) ) - S; 

% loop 
err = 1e-6; 
while abs(lat - lat0) > 1e-4

    C = R_E / sqrt( 1 - e^2 * ( sin(lat) )^2 ); 
    S = R_E * (1-e^2) / sqrt( 1 - e^2*( sin(lat) )^2 ); 

    tan_lat = (r_z + C * e^2 * sin(delta)) / r_delta; 
    lat = atan(tan_lat); 
    if abs(90 - i) > 1
        h = r_z / ( cos(delta) ) - C; 
    else
        h = r_z / ( sin(delta) ) - S; 
    end 

end 

% h = h + R_E; 
lla = [lat; lon; h]; 

end 



% ------------------------------------------------------------------------
% ECEF2LLA - convert earth-centered earth-fixed (ECEF)
%            cartesian coordinates to latitude, longitude,
%            and altitude
%
% USAGE:
% [lat,lon,alt] = ecef2lla(x,y,z)
%
% lat = geodetic latitude (radians)
% lon = longitude (radians)
% alt = height above WGS84 ellipsoid (m)
% x = ECEF X-coordinate (m)
% y = ECEF Y-coordinate (m)
% z = ECEF Z-coordinate (m)
%
% Notes: (1) This function assumes the WGS84 model.
%        (2) Latitude is customary geodetic (not geocentric).
%        (3) Inputs may be scalars, vectors, or matrices of the same
%            size and shape. Outputs will have that same size and shape.
%        (4) Tested but no warranty; use at your own risk.
%        (5) Michael Kleder, April 2006
function [lla] = ecef2lla_2(ecef)

x = ecef(1); y = ecef(2); z = ecef(3); 

% WGS84 ellipsoid constants:
a = 6378137;
e = 8.1819190842622e-2;
% calculations:
b   = sqrt(a^2*(1-e^2));
ep  = sqrt((a^2-b^2)/b^2);
p   = sqrt(x.^2+y.^2);
th  = atan2(a*z,b*p);
lon = atan2(y,x);
lat = atan2((z+ep^2.*b.*sin(th).^3),(p-e^2.*a.*cos(th).^3));
N   = a./sqrt(1-e^2.*sin(lat).^2);
alt = p./cos(lat)-N;
% return lon in range [0,2*pi)
lon = mod(lon,2*pi);
% correct for numerical instability in altitude near exact poles:
% (after this correction, error is about 2 millimeters, which is about
% the same as the numerical precision of the overall function)
k=abs(x)<1 & abs(y)<1;
alt(k) = abs(z(k))-b;

lla = [lat; lon; alt]; 

end 


