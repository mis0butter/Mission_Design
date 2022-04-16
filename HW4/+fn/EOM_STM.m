function dXSTM = EOM_STM(et, XSTM, Amat_fn, n)
% ------------------------------------------------------------------------
% Purpose: Generate EOM for satellite orbiting earth due to gravity (EGM 96),
% lunisolar perturbations, SRP, and drag 
% 
% Inputs 
%   t  = [7x1] time (ET epoch) vector 
%   X  = [49x1] state vector in ECI frame (inertial) 
% 
% Outputs 
%   dX = [49x1] derivative of state vector 
% ------------------------------------------------------------------------

% force column vector 
dXSTM = zeros(length(XSTM), 1);   

% EOM 
dX         = fn.EOM(et, XSTM); 
dXSTM(1:n) = dX(1:n); 

% STM stuff 
if n == 6
    Amat = Amat_fn(XSTM(1), XSTM(2), XSTM(3), XSTM(4), XSTM(5), XSTM(6)); 
else
    Amat = Amat_fn(XSTM(1), XSTM(2), XSTM(3), XSTM(4), XSTM(5), XSTM(6), XSTM(7)); 
end
STM  = XSTM(n+1:end); 
STM  = reshape(STM, [n n]); 
dSTM = Amat*STM; 
dSTM = reshape(dSTM, [n^2, 1]); 

dXSTM(n+1:end) = dSTM; 

end 



