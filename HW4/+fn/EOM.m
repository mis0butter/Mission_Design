function dX = EOM(et, X)
% ------------------------------------------------------------------------
% Purpose: Generate EOM for satellite orbiting earth due to geopotential,
% lunisolar, SRP, and drag perturbations 
% 
% Inputs 
%   t   = [1x1] time (ET epoch) vector 
%   X   = [6x1] state vector in ECI frame (inertial) 
% 
% Outputs 
%   dX  = [6x1] derivative of state vector 
% ------------------------------------------------------------------------

% global const 
% sun mu (m^3/s^2)
mu_sun_m = 1.32712440018e20; 
mu_sun_km = mu_sun_m / (1000^3); 
mu = mu_sun_km; 

% force column vector. Check if X is numeric or symbolic 
n = length(X); 
dX = zeros(n, 1);   

% ------------------------------------------------------------------------
% Set velocity and CD 
dX(1:3) = X(4:6);

% ------------------------------------------------------------------------
% % accel due to gravity
% g = fn.g_J2J3J4(X); 
% dX(4:6) = dX(4:6) + g; 

% accel due to point mass (not needed when geopotential gravity is present)
r       = norm(X(1:3)); 
r_norm  = sqrt( X(1)^2 + X(2)^2 + X(3)^2 ); 
dX(4:6) = ( - mu / r_norm^3 ) * X(1:3); 

% ------------------------------------------------------------------------
% % accel due to lunisolar perturbation 
% [a_sun, a_moon] = fn.lunisolar(et, X); 
% dX(4:6) = dX(4:6) + a_sun + a_moon; 

% ------------------------------------------------------------------------
% % accel due to SRP 
% a_srp   = fn.a_SRP(et, X); 
% dX(4:6) = dX(4:6) + a_srp; 

% ------------------------------------------------------------------------
% % CHECK UNITS. USE KM !!!
% % accel due to drag 
% pA      = const.p0 * exp( -(r - const.r0_drag)/const.H ); 
% VA      = X(4:6) - cross( [0; 0; const.wN], X(1:3) ); 
% VAnorm  = norm(VA); 
% a_drag  = - 1/2 * const.CD * const.A / const.m_SC * pA * VAnorm * VA;
% a_drag  = a_drag / 26.5;  % Correct to match Jah's Amat??? 
% dX(4:6) = dX(4:6) + a_drag; 

end 