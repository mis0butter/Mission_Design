function drv = TwoBod_6states(t, rv)
% ------------------------------------------------------------------------
% Inputs 
%   t = [Nx1] time vector (orbit is Keplerian, doesn't matter) 
%   x = [6x1] state vector 
% 
% Outputs 
%   dx = [6x1] derivative of state vector 
% ------------------------------------------------------------------------

global mu 

drv = zeros(6, 1);   % force column vector 

% dx1 = x4 
% dx2 = x5 
% dx3 = x6 
% dx4 = (-u/r^3) * x1
% dx5 = (-u/r^3) * x2 
% dx6 = (-u/r^3) * x3 

drv(1:3) = rv(4:6); 
r_norm  = sqrt( rv(1)^2 + rv(2)^2 + rv(3)^2 ); 
drv(4:6) = ( - mu / r_norm^3 ) * rv(1:3); 

end 