function [rv] = oe2rv(oe)
% ------------------------------------------------------------------------ 
% Purpose: Convert orbital elements and time past epoch to the classic 
% Cartesian position and velocity
% 
% Inputs: 
%   oe      = [6x1] or [1x6] orbital elements 
%   delta_t = t - t0 time interval 
%   mu      = Gravity * Mass (of Earth) constant 
% 
% Outputs: 
%   rv      = position and velocity state vector 
% ------------------------------------------------------------------------ 

% global delta_t 
global mu 

a       = oe(1); 
e       = oe(2); 
i       = oe(3); 
w       = oe(4); 
O     = oe(5); 
% M0      = oe(6); 
nu      = oe(6); 

% nu is TRUE ANOMALY --> use Kepler's to calculate MEAN ANOMALY 
% E = 2*atan( sqrt( (1-e)/(1+e) ) * tan(nu/2) ); 
% M = M0 + sqrt( mu/a^3 ) * (delta_t); 
% E = keplerEq(M, e, eps); 
% E = kepler(M, e); 
% nu = 2*atan( sqrt( (1+e)/(1-e) ) * tan(E/2) ); 

p = a * ( 1 - e^2 );            % intermediate variable 
r = p / ( 1 + e*cos(nu) );      % r_magnitude, polar coordinates 

% Perifocal position and velocity 

r_pf = [ r * cos(nu); r * sin(nu); 0 ]; 
v_pf = [ -sqrt(mu/p) * sin(nu); sqrt(mu/p) * (e + cos(nu)); 0 ]; 

% Perifocal to ECI transformation, 3-1-3 rotation 
R11 = cos(O)*cos(w) - sin(O)*sin(w)*cos(i); 
R12 = -cos(O)*sin(w) - sin(O)*cos(w)*cos(i); 
R13 = sin(O)*sin(i); 

R21 = sin(O)*cos(w) + cos(O)*sin(w)*cos(i); 
R22 = -sin(O)*sin(w) + cos(O)*cos(w)*cos(i); 
R23 = -cos(O)*sin(i); 

R31 = sin(w)*sin(i); 
R32 = cos(w)*sin(i); 
R33 = cos(i); 

R = [R11 R12 R13; R21 R22 R23; R31 R32 R33]; 

% Transform perifocal to ECI frame 
r_vec = R * r_pf; 
v_vec = R * v_pf; 

% Position and state vector 
rv = [r_vec; v_vec]; 

end 

%% Kepler equation solvers 

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

function E = kepler(M, e)
    f = @(E) E - e * sin(E) - M;
    E = fzero(f, M);  % <-- I would use M as the initial guess instead of 0
end