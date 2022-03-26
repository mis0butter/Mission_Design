% ASE 389 Orbit Determination
% HW 2
% Junette Hsin 

clear; 

% positionA = [100 100 600 600]; 
% positionB = [100 100 700 900]; 

%% Problem 1a: 
% Assume an orbit plane coordinate system with a gravitational parameter of 1, i.e., µ = 1. The
% equations of motion are:
% ddx = - x / r^3
% ddy = - y / r^3 
% r^2 = x^2 + y^2 
% 
% Generate a “true” solution by numerically integrating the equations of motion for the initial
% conditions:
% X(t0) = [x; y; dx; dy] = [1; 0; 0; 1] 
% Save the values of the state vector X(ti) for ti = i · 10 time units (TU); i = 0, . . . , 10.
% Provide X(ti) for t1 and t10 in the writeup.
% In your write-up, please indicate which integrator you used, what the tolerance was set
% to, and any other details necessary. Note, if you use a fixed time-step integrator, set the
% time-step to be smaller than 10 TU, but only save the data at 10 TU intervals.

global mu 
mu = 1; 

% set ode45 params 
rel_tol = 3e-14;         % 1e-14 accurate; 1e-6 coarse 
abs_tol = 1e-16; 
options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

rv0 = [1; 0; 0; 1]; 
dt = 0.01; 

% integrate 
[t, rv] = ode45(@fn.TwoBod_4states, [0:dt:100], [rv0], options); 

%% Problem 1b: 
% Perturb the previous set of initial conditions by an amount
% X∗(t0) = X(t0) − δX(t0)
% (notice that the perturbation is subtracted!), where 
% δX(t0) = [1e-6; -1e-6; 1e-6; 1e-6] 

drv0 = [1e-6; -1e-6; 1e-6; 1e-6]; 
STM0 = eye(4); 
STM0 = reshape(STM0, [16 1]); 

rvSTM0 = [rv0 - drv0; STM0]; 

% integrate 
[tstar, rvstar] = ode45(@fn.TwoBod_4states_STM, [0:dt:100], [rvSTM0], options); 

STMf = rvstar(end, 5:20); 
STMf = reshape(STMf, [4 4]); 

%% Problem 1c: 
% For this problem, Φ(ti, t0) is symplectic. Demonstrate this for Φ(t10, t0) by multiplying it by
% Φ^−1(t10, t0), given by Eq. 4.2.22 in the text. Provide Φ^−1(t10, t0) and show that the product
% with Φ(t10, t0) is the identity matrix.

STMf1 = STMf(1:2, 1:2); 
STMf2 = STMf(1:2, 3:4); 
STMf3 = STMf(3:4, 1:2); 
STMf4 = STMf(3:4, 3:4); 

STMfinv = [ STMf4',  -STMf2'; ... 
             -STMf3', STMf1' ]; 
         
STMfinv * STMf 

%% Problem 1d: 
% Calculate the perturbation vector, δX(ti), by the following methods:
% (1) δX(ti) = X(ti) − X∗(ti)
% (2) δX(ti) = Φ(ti, t0)δX(t0)
% and compare the results of (1) and (2). Provide the numeric results of (1) and (2) at t1 and
% t10 in the write-up, along with δX(ti) − Φ(ti, t0)δX(t0). How closely do they compare?

% t1 = 10 TU
i = 10 / dt + 1; 
STMi = rvstar(i, 5:20); 
STMi = reshape(STMi, [4 4]); 

drv1 = rv(i,:) - rvstar(i,1:4); 
drv2 = STMi * drv0; 

ddrvt1 = drv1' - drv2; 

% t10 = 100 TU
i = 100 / dt + 1; 
STMi = rvstar(i, 5:20); 
STMi = reshape(STMi, [4 4]); 

drv1 = rv(i,:) - rvstar(i,1:4); 
drv2 = STMi * drv0; 

ddrvt10 = drv1' - drv2; 

%% Problem 2: Given the observation state relation y = H x + eps, where x is a scalar and
% y = [1; 2; 1]
% W = [2 0 0; 0 1 0; 0 0 1]; 
% H = [1; 1; 1]; 
% with a priori information xbar = 2 and Wbar = 2:
% 
% Problem 2a: Using the batch processing algorithm, what is xˆ? In the write-up, outline the method
% employed in the code.

% W matrix is inv(R) 
% Wbar = inv(P) 

% observation states 
y = [1; 2; 1]; 
W = [2 0 0; 0 1 0; 0 0 1]; 
H = [1; 1; 1]; 

x0 = 2; 
W0 = 2; 

% Lambda = inv(P) = W0  
P0      = inv(W0); 
Lambda  = W0; 
R       = inv(W); 

% N = inv(P) * x0 = W0 * x0 
N = W0 * x0; 

% accumulate 
Lambda  = Lambda + H' * W * H; 
N       = N + H' * W * y; 

% normal equation 
xhat = inv(Lambda) * N; 

% xhat = inv(H'*inv(R)*H + inv(P)) * H'*inv(R)*y + inv(P)*x0; 
% xhat = inv(H' * W * H + W0) * H' * W + W0*x0; 

%%

x0 = 2;
W0 = 2; % Pinv
W = [2 0 0;0 1 0;0 0 1]; % Rinv
H = [1 1 1]';

lambda = W0;
N = W0*x0;
y = [1 2 1]';
lambda = lambda + H'*W*H;
N = N + H'*W*y;
xhat = N/lambda;
eps = y - H*xhat;

%% Problem 2b: What is the best estimate of the observation error, eps? 

e = y - H*xhat; 


%%     
%% subfunctions 


