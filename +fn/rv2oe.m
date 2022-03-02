function [oe, Tp] = rv2oe(rv)
% ------------------------------------------------------------------------
% Inputs 
%   rv = [6x1] position and velocity states vector in ECI frame 
% 
% Outputs 
%   oe = [6x1] orbital elements: a, e, i, w, Omega, nu
%           a   = semimajor axis 
%           e   = eccentricity 
%           i   = inclination 
%           w   = argument of perigee 
%           O   = right ascension of ascending node 
%           nu  = true anomaly 
%           Tp  = time of perigee passing  (optional output) 
% ------------------------------------------------------------------------

global mu 

r = rv(1:3); 
v = rv(4:6); 

% angular momentum 
h       = cross(r,v); 

% node vector 
nhat    = cross( [0 0 1], h ); 

% eccentricity 
evec    = ( (norm(v)^2 - mu/norm(r))*r - dot(r,v)*v ) / mu; 
e       = norm(evec); 

% specific mechanical energy 
energy  = norm(v)^2/2 - mu/norm(r); 

% semi-major axis and p
if abs(e-1.0)>eps
   a = -mu/(2*energy); 
   p = a*(1-e^2); 
else
   p = norm(h)^2/mu; 
   a = inf; 
end

% inclination 
i = acos(h(3)/norm(h)); 

% right ascension of ascending node (check for equatorial orbit) 
if i > 0.000001
    O = acos( nhat(1)/norm(nhat) ); 
else
    O = 0; 
end
if isnan(O)
    O = 0; 
end
if nhat(2)<0
   O = 2*pi - O; 
end

% argument of perigee 
if e > 0.000001
    w = acos(dot(nhat,evec)/(norm(nhat)*e)); 
else
    w = 0; 
end
if isnan(w)
    w = 0; 
end
% if e(3)<0
%    argp = 360-argp
% end

% true anomaly 
nu = acos( dot(evec,r) / (e*norm(r)) );  
% if dot(r,v)<0
%    nu = 360 - nu
% end

%Apply Quadrant Checks to All Determined Angles
% idx = nhat(2) < 0; if any(idx); O(idx) = 2*pi - O(idx);  end
% idx = evec(3) < 0; if any(idx); w(idx) = 2*pi - w(idx); end
idx = dot(r,v) < 0; if any(idx); nu(idx) = 2*pi - nu(idx); end

% RAAN = atan2d(h(1),-h(2));
% 
% ArgLat --> arg of perigee
% ArgLat = atan2d(R0(3)/sind(i),R0(1)*cosd(RAAN) + R0(2)*sind(RAAN));
% 
% TA = atan2d(sqrt((a*(1-e^2))/(mu))*dot(V0,R0),(a*(1-e^2))-norm(R0));

oe = [a; e; i; w; O; nu]; 

% perigee passing time 

rnorm = norm(r); 
n = sqrt(mu/a^3); 
E = acos( rnorm/a * cos(nu) + e );
M = E - e*sin(E); 
Tp = M/n; 

end