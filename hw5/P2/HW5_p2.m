
clear; clc;
format long g

%{
Some things:
You'll need OE2RV and RV2OE
It's currently set up for
    OE = [a, e, i, LAN, omega, nu]
I started with EGM96, and was going to change that to EIGEN.. whatever, at a later time..
not sure if it matters.
    EGM96_MOD is a slight variation to EGM96 that Jah used
    I load in EGM96 once instead of needing to load it in
    each function call for gravitysphericalharmonic

Section 3.2 should work using "run_sso"

The idea in Section 3.3 is to vary initial e such that the 
ex and ey components oscillate around a single value 
instead of precess.
    I have not included any automatic variation of e

I don't know to get a proper omega, LAN, and M0 using Eq. 11.
    Eq.11 is implemented via function "iterate" in this code
    gradient descent


%}


% path for mex'd gravitysphericalharmonic
addpath 'C:\Users\User\Documents\Spectre\2022\Spring\MissionDesign\HW\HW5\codegen\mex\gravitysphericalharmonic_fast'
global EGM96_MOD degree
EGM96_MOD = load('EGM96_MOD.mat');
% warning('off','last');


mu = 398600.4415; % km^3/s^2
w = 7.292115146706979e-05;

global JD0 lon_N tol_prop
JD0 = 2453832.168203696; % 2006, 4, 6.668203696142882
lon_N = 28.919887923739253;

% Semimajor Axis
% Eccentricity
% Inclination
% Right Ascention of the Ascending Node
% Argument of Periapse
% True Anomaly

% D'Amico conditions at SGP4 solution (Table 1):
a = 6892.950; %6883.509411974535;
e = 0;
inc = 97.4464*pi/180; % 1.7003362447800794;
LAN = 104.2749*pi/180;
omega = 180.0270*pi/180;
M0 = 180.0965*pi/180;
nu0 = M0;
OE_cycle0 = [a, e, inc, LAN, omega, nu0]';

% Tolerance/degree for simulation accuracy
tol_prop = 1e-9;
degree = 0;

% Options to fill out Table 2 in D'Amico
run_sso = 1; % sun-synchronous
run_frozen = 1;

if 0
    % runtime testing
    % tic
    degree = 20;
    tol_prop = 1e-9; %2.22045e-14;
    [OE, err_lon, err_lat] = prop_cycle(OE_cycle0)
    toc
end

if run_sso
    % sun-sync condition
    %{
    iterate
        function to implement Gradient Descent
        takes in OE_cycle0, the current, initial
        OEs (a, e, i, LAN, omega, nu), max iterations,
        a scale factor on the rate of descent (alpha),
        and the region over which the functions f and g
        are considered "linear", da, and di
    
    Could try to do this in 1 large iteration? There might
    be faster ways to make it work.
    As of now, I have an iteration for the iteration..
        i = 1:num_iterations
        every i, lower alpha (rate of descent), and
        the linear region defined by da and di
        algorithm -does- converge slowly to highly accurate a and inclination

    save()
        save a copy of the OEs, so you don't have to run this
        every time (change run_sso and run_frozen flags)
    %}

    num_iterations = 2;
    max_iter = 2;
    alpha = 1.0;
    da = 0.01;
    di = 1e-7;
    for i = 1:num_iterations
        fprintf('iteration %d\n', i);
        OE_cycle0 = iterate(OE_cycle0, max_iter, alpha, da, di);
        da = da/10;
        di = di/10;
        alpha = alpha/2;
    end

%     num_iterations = 10;
%     max_iter = 3;
%     alpha = 1.0;
%     da = 1.0;
%     di = 1e-5;
%     for i = 1:num_iterations
%         fprintf('iteration %d\n', i);
%         OE_cycle0 = iterate(OE_cycle0, max_iter, alpha, da, di);
%         da = da/5;
%         di = di/5;
%         alpha = alpha/2;
%     end

    save('OE_cycle0.mat', 'OE_cycle0');
end

if run_frozen
    % frozen orbit condition
    %{
    The idea is to propagate forward 11-day repeats, 10 times (110 days),
    and output eccentricity each time. Then, choose e~=0 as the initial
    eccentricity, and try again.. until the propagated e becomes "small"
    (see Section 3.2 D'Amico)

    Either continue with SSO OE's or load them up
    Assumes tolerance and degree are the same as in SSO

    num_cycles
        the paper loops over 10 cycles, or 10 11-day propagations
        - using the initial conditions from SSO, propagate 110 days
        - save OEs after 11 days, and restart propagator
    
    err_lon, err_lat (rad, but printed ind deg)
        the error in the 11-day repeat
        ideally this is zero after every cycle
    
    
    %}
    if ~run_sso
        data = load('OE_cycle0.mat');
        OE_cycle0 = data.OE_cycle0;
    end
%     OE_cycle0 = [6893.03379670928;
%                          0;
%           1.70037842759382;
%           1.81994033218783;
%           3.14206389248783;
%           3.14327689631797];
%     %

    num_cycles = 10;
    OE_cycle = zeros(6,num_cycles);
    for k = 1:num_cycles
        fprintf('cycle %d\n', k);
        [OE, err_lon, err_lat] = prop_cycle(OE_cycle0); % propagate 11 days
        disp([OE(1), err_lon*180/pi, err_lat*180/pi]); % display cycle error
        OE_cycle0 = OE; % update OEs
        OE_cycle(:,k) = OE; % save OEs
        JD0 = JD0 + 11; % update initial propagation time
    end

    % build ex and ey def. in Section 3.2 of D'Amico
    e_vec = OE_cycle(2,:)';
    omega_vec = OE_cycle(5,:)';
    ex = e_vec .* cos(omega_vec);
    ey = e_vec .* sin(omega_vec);

    figure;
    plot(ex, ey);
    xlabel('ex');
    ylabel('ey');

end


function OE_out = iterate(OE_in, max_iter, alpha, da, di)

% Gradient descent
% Determines Jacobian from the prior state, updates state
%   There's probably a faster way to do this
%   It requires 5 propagations over 11 days to determine J and dz :/

tol = 1e-6;
i = 0;

x = [OE_in(1); OE_in(3)];
while(i < max_iter)

    % derivative w.r.t. a
    [~, df1, dg1] = prop_cycle([x(1)-da/2; OE_in(2); x(2); OE_in(4:end)]);
    [~, df2, dg2] = prop_cycle([x(1)+da/2; OE_in(2); x(2); OE_in(4:end)]);
    dfda = (df2-df1)/da;
    dgda = (dg2-dg1)/da;

    % derivative w.r.t. inc
    [~, df1, dg1] = prop_cycle([x(1); OE_in(2); x(2)-di/2; OE_in(4:end)]);
    [~, df2, dg2] = prop_cycle([x(1); OE_in(2); x(2)+di/2; OE_in(4:end)]);
    dfdi = (df2-df1)/di;
    dgdi = (dg2-dg1)/di;

    J = [dfda dfdi; dgda dgdi]; % Jacobian

    % determine current error, solve for dlam, dphi
    [~, err_lon, err_lat] = prop_cycle([x(1); OE_in(2); x(2); OE_in(4:end)]);
    dz = [err_lon; err_lat];
    dx = -inv(J) * dz;

    dx0 = alpha*dx; % multiply by alpha (rate of descent) if necessary
    disp([i dx0(1) dx0(2)]);

    x = x + dx0; % update state
    eps = norm(dx0);
    if(eps < tol)
        break;
    end
    
    i = i + 1;
end

OE_out = zeros(6,1);
for j = 1:6
    if j == 1
        OE_out(j) = x(1);
    elseif j == 3
        OE_out(j) = x(2);
    else
        OE_out(j) = OE_in(j);
    end
end

end


function R_ECI_ECF = eci_to_ecf_ss(JD)
% simple spinner Earth
D = JD - 2451545.0;
theta_g = 18.697374458 + 24.06570982*D;
GMST2 = wrapTo180(theta_g*15);
R_ECI_ECF = get_Rz(GMST2*pi/180);
end

function Rz = get_Rz(theta)
% passive rotation z
R1 = [cos(theta) -sin(theta) 0];
R2 = [sin(theta) cos(theta) 0];
R3 = [0 0 1];
Rz = [R1; R2; R3]';
end

function M0 = nu2M0(nu, a, e, t, t0, mu)
E = 2*atan(tan(nu/2)/(sqrt((1+e)/(1-e))));
M = E - e*sin(E); % since nu is at t0
n = sqrt(mu/a^3);
M0 = M - n*(t-t0);
end

function nu = M02nu(a, e, M0, t, t0, mu)
n = sqrt(mu/a^3);
M = M0 + n*(t-t0);
E = calc_E(e, M);
nu = 2*atan(sqrt((1+e)/(1-e)) * tan(E/2.0)); % rad
end


function [OE, err_lon, err_lat] = prop_cycle(OE_cycle0)

% propagate forward 11 days from JD0 (global var)
% output propagated OE, and errors in lon and lat at 11-day repeat

global JD0 lon_N tol_prop
mu = 398600.4415; % km^3/s^2

dT = 11 * 86400;

% tol = 2.22045e-14;
tol = tol_prop;
options = odeset('RelTol', tol, 'AbsTol', tol);

% [a, e, inc, LAN, omega, M0]
% nu = M02nu(a, e, M0, 0, 0, mu);

y0 = COEstoRV(OE_cycle0, mu);
[~, y] = ode45(@ode_EKF, [0 dT], y0, options); % , JD0, table, N_data, degree);
yf = y(end,:)';
r = yf(1:3);
v = yf(4:6);

% OE = cm_ut.RV2OE(r, v, 0, 0, mu)
% LAN_prop = OE[4]
OE = RVtoCOEs(r,v,mu);
LAN_prop = OE(4);

lat = asin(r(3)/norm(r))*180/pi;
JD = JD0 + dT/86400;
D = JD - 2451545.0;
theta_g = 18.697374458 + 24.06570982*D;
GMST2 = wrapTo180(theta_g*15);
lon_N2 = LAN_prop*180/pi - GMST2;

err_lon = abs(lon_N2-lon_N) * pi/180;
err_lat = abs(lat-0) * pi/180;

end



function y_dot = ode_EKF(tau, y)

global EGM96_MOD JD0 degree
% 20x20
% luni-solar
% total s/c SRP
% cylindrical shadow

% if need to change, run EOM() with new potential values
% copy and paste v_dot_s into v_dot
%   careful that r_mag2 is used, not r_mag

r = y(1:3);
v = y(4:6);

rx = r(1);
ry = r(2);
rz = r(3);
% vx = v(1);
% vy = v(2);
% vz = v(3);

% gravitational
mu = 398600.4415; % km^3/s^2
R = 6378.1363; % earth radius, km
J2 = 0.00108248; %0.0010826267;
J3 = -0.0000025327;

r_mag = norm(r);
r_mag2 = r_mag^2;
r_dot = v;

if degree == 0
%     v_dot_g = -mu/r_mag^3 * r; % no gravitational perturbation
    v_dot_g = [ -(15*J3*mu*rx*R^3*r_mag2*rz - 35*J3*mu*rx*R^3*rz^3 + 3*J2*mu*rx*R^2*r_mag2^2 - 15*J2*mu*rx*R^2*r_mag2*rz^2 + 2*mu*rx*r_mag2^3)/(2*r_mag2^(9/2)), ... 
            -(15*J3*mu*ry*R^3*r_mag2*rz - 35*J3*mu*ry*R^3*rz^3 + 3*J2*mu*ry*R^2*r_mag2^2 - 15*J2*mu*ry*R^2*r_mag2*rz^2 + 2*mu*ry*r_mag2^3)/(2*r_mag2^(9/2)), ...
            -(- 3*J3*mu*R^3*r_mag2^2 + 30*J3*mu*R^3*r_mag2*rz^2 - 35*J3*mu*R^3*rz^4 + 9*J2*mu*R^2*r_mag2^2*rz - 15*J2*mu*R^2*r_mag2*rz^3 + 2*mu*r_mag2^3*rz)/(2*r_mag2^(9/2))]';
    %

else
    JD1 = JD0 + tau/86400;
    % R_ECF_ECI = eci_to_ecf(JD1, table, N_data);
    R_ECF_ECI = eci_to_ecf_ss(JD1);
    r_ECF = R_ECF_ECI * r;

    % gravitysphericalharmonic_fast_mex
    %   made faster than normal gravitysphericalharmonic b/c
    %   this takes in the EGM96 constants as variables, rather
    %   than loading them on each function call

    [gx, gy, gz] = gravitysphericalharmonic_fast_mex(r_ECF'*1e3, EGM96_MOD.GM, EGM96_MOD.Re, degree, EGM96_MOD.C, EGM96_MOD.S);
    v_dot_ECF = [gx gy gz]'/1e3;
    
    v_dot_g_new = R_ECF_ECI' * v_dot_ECF;

    v_dot_g = v_dot_g_new;
end

v_dot = v_dot_g;

y_dot = [r_dot; v_dot];

end

