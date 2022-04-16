function [a_sun, a_moon] = lunisolar(et, X)
% From Born 
% Perturbed Motion Equation 2.3.39 
% modified for Neptune-Triton system 

global const 

% Moon wrt Earth 
target      = 'Moon';
frame       = 'J2000';
observer    = 'Earth';
abcorr      = 'NONE';

% Triton wrt Neptune 
target      = 'Triton'; 
frame       = 'J2000'; 
observer    = 'Neptune'; 
abcorr      = 'NONE'; 

% get states --> Neptune to Triton 
X_Nep_Tri = fn.spice_state(et, target, frame, abcorr, observer); 
X_Nep_Tri = X_Nep_Tri'; 

% Sun wrt Neptune 
target  = 'Sun'; 

% get states --> Neptune to Sun 
X_Esun = fn.spice_state(et, target, frame, abcorr, observer); 
X_Esun = X_Esun'; 

% pos vector --> sat to sun 
X_sunsat = X(1:6) - X_Esun;
X_satsun = -X_sunsat; 

% Relative position vector of satellite wrt Triton 
X_moonsat = X(1:6) - X_Nep_Tri; 
X_satmoon = -X_moonsat; 

% lunisolar acceleration 
a_sun  = const.muS * ( X_satsun(1:3)/norm(X_satsun(1:3))^3 - X_Esun(1:3)/norm(X_Esun(1:3))^3 );  
a_moon = const.muT * ( X_satmoon(1:3)/norm(X_satmoon(1:3))^3 - X_Nep_Tri(1:3)/norm(X_Nep_Tri(1:3))^3 ); 

end



