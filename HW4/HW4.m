%% Junette Hsin 

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
O0 = 0 * pi/180;         % deg --> rad
M0 = 39.7178 * pi/180;   % deg --> rad (should be true anomaly) 

% mean motion --> semimajor axis 
n = 15.50094 / 86400 * (2*pi);   % rev/day --> rad/s 
mu_E_m3 = 3.986004418e14;  % m^3/s^2
mu_E_km3 = mu_E_m3 * ( 1e-3 )^3 ; 
a0 = (mu_E_km3 / n^2)^(1/3); 

oe0 = [a0; e0; i0; w0; O0; M0]; 
rv0 = rvOrb.orb2rv(oe0, mu_E_km3); 

% set ode45 params 
rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-10; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

% propagate orbit 
[t, rv] = ode45(@fn.EOM, [0 : 15*60], rv0, options); 

lla = []; 
lla_1 = []; 
for i = 1:length(rv)
    
    lla(i,:) = ecef2lla(rv(i,1:3)); 
    lla_1(i,:) = ecef2lla_1(rv(i,:)', R_E, mu_E_km3); 
    lla_1(i, 1:2) = lla_1(i, 1:2) * 180/pi; 
    
end 
% fn.plot3_xyz(rv); 

% figure()
%     plot(lla(:,2), lla(:,1), 'b')
%     hold on; grid on; 
%     plot(lla(1,2), lla(1,1), 'bo')
%     plot(lla(end,2), lla(end,1), 'b^')
%     
%     plot(lla_1(:,2), lla_1(:,1), 'r')
%     plot(lla_1(1,2), lla_1(1,1), 'ro')
%     plot(lla_1(end,2), lla_1(end,1), 'r^')
%     
%     xlabel('longitude'); ylabel('latitude'); 

figure()
lh = plot(lla_1(:, 2), lla_1(:, 1), 'r'); 
hold on; grid on; 
xl = xlim; yl = ylim; 

delete(lh); 
xlim(xl); ylim(yl); 
for i = 1:length(lla_1)

    cla 

%     % mathworks ecef2lla 
%     plot(lla(1,2), lla(1,1), 'bo')
%     plot(lla(1:i,2), lla(1:i,1), 'b', 'linewidth', 1.5)
%     plot(lla(end,2), lla(end,1), 'b^')
    
    % my ecef2lla 
    plot(lla_1(1,2), lla_1(1,1), 'ro'); 
    plot(lla_1(1:i, 2), lla_1(1:i, 1), 'r--'); 
    lh = plot(lla_1(i, 2), lla_1(i, 1), 'r^'); 
    
    title(sprintf( '%d / %d', i, length(lla_1) )); 
    pause(0.001) 
    
end 
% plot(lla_1(end,2), lla_1(end,1), 'r^')

%% test 

rv_test = [6524.834; 6862.875; 6448.296]; 
lla_1 = ecef2lla_1(rv_test, oe0, R_E); 
lla_3 = ecef2lla_3(rv_test, oe0, R_E); 
lla = ecef2lla(rv_test')';
lla(1:2) = lla(1:2) * pi/180; 
lla_2 = ecef2lla_2(rv_test); 



%% subfunctions 

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



%% 

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


