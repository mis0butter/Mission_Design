% ASE 389 Orbit Determination
% HW 1
% Junette Hsin 

%% Problem 1 

global mu 

mu = 398600.5 ; 
r  = [ -2436.45; -2436.45; 6891.037 ]; 
v  = [ 5.088611; 5.088611; 0 ]; 

rv = [r; v]; 
oe = rv2oe(rv);  

%% Problem 2 

rv = oe2rv(oe); 

%% Problem 3 

x = -2436.45; 
y = -2436.45; 
z = 6891.037; 
mu = 398600.5; 

dux = -2*mu*x / ( x^2 + y^2 + z^2 )^2; 
duy = -2*mu*y / ( x^2 + y^2 + z^2 )^2;
duz = -2*mu*z / ( x^2 + y^2 + z^2 )^2;

rnorm = sqrt( x^2 + y^2 + z^2 ); 

dux = -mu*x / ( rnorm )^3; 
duy = -mu*y / ( rnorm )^3; 
duz = -mu*z / ( rnorm )^3;

%% Problem 4 

a = oe(1); 
T   = abs(2 * pi * sqrt(a^3 / mu));        % period 

toler   = 1e-8;         % 1e-14 accurate; 1e-6 coarse 
options = odeset('reltol', toler, 'abstol', toler ); 
[t,x] = ode45(@TwoBod_6states, [0 2*T], [r; v], options); 

for i = 1:length(t)
    rnorm(i) = norm(x(i, 1:3)); 
    vnorm(i) = norm(x(i, 4:6)); 
    H(i, :) = cross(x(i, 1:3), x(i, 4:6)); 
    hnorm(i) = norm(H(i, :)); 
end 

anorm = 0; 
for i = 2:length(t)
    a = (x(i, 4:6) - x(i-1, 4:6)) / ( t(i) - t(i-1) ); 
    anorm(i) = norm(a); 
end 

% ------------------------------------------------------------------------

name = 'Problem 4: 2-Body EOM'; 
h = figure('name', name); 

    % position 
    subplot(3,1,1)
        plot(t, rnorm); grid on 
        title('r norm') 
        ylabel('km')

    % velocity 
    subplot(3,1,2) 
        plot(t, vnorm); grid on 
        title('v norm') 
        ylabel('km/s')

    % acceleration 
    subplot(3,1,3) 
        plot(t, anorm); grid on 
        title('a norm'); 
        ylabel('km/s^2')
        xlabel('time (sec)') 

        sgtitle(name)

save_pdf(h, 'prob4_2bodeom'); 

% ------------------------------------------------------------------------

name = 'Problem 4: 2-Body EOM Orbit'; 
h = figure('name', name); 
    plot3(x(:,1), x(:,2), x(:,3)); hold on; grid on; 
    plot3(x(1,1), x(1,2), x(1,3), 'o')
    plot3(x(end,1), x(end,2), x(end,3), 'x') 
    xlabel('x (km)')
    ylabel('y (km)') 
    zlabel('z (km)') 
    legend('orbit', 'start', 'end')
    
    sgtitle(name) 

save_pdf(h, 'prob4_2bodeom_orbit'); 

% ------------------------------------------------------------------------

clear oe 
for i = 1:length(t)
    oe(i,:) = rv2oe(x(i,:)); 
end 

labels = {'a', 'e', 'i', '\omega', '\Omega', '\nu'}; 
units = {'km', '', 'rad', 'rad', 'rad', 'rad'}; 
name = 'Problem 4: 2-Body Orbital Elements'; 
h = figure('name', name, 'position', [100 100 500 600]); 
    for i = 1:6
        subplot(6,1,i)
        plot(t, oe(:, i)); grid on 
        title(labels{i}); 
        ylabel(units{i}); 
    end 
    xlabel('time (sec)') 
    sgtitle(name)

save_pdf(h, 'prob4_2bodoes'); 

% ------------------------------------------------------------------------

name = 'Problem 4: 2-Body Specific Angular Momentum'; 
h = figure('name', name); 
subplot(2,1,1) 
    scatter3(H(:,1), H(:,2), H(:,3)); grid on 
    xlabel('x (km^2/s)')
    ylabel('y (km^2/s)')
    zlabel('z (km^2/s)') 
    title('h (scatter plot)') 

subplot(2,1,2) 
    plot(t, hnorm); grid on 
    xlabel('time (sec)') 
    ylabel('km^2/s') 
    title('h norm vs time')
    
    sgtitle(name) 
    
save_pdf(h, 'prob4_angmom')

%% Problem 5 

% specific kinetic energy 
for i = 1:length(t) 
    T(i) = 0.5 * vnorm(i)^2; 
    U(i) = mu / rnorm(i); 
end 
E = T - U; 

name = 'Problem 4: 2-Body Specific Energy'; 
h = figure('name', name); 
    subplot(2,1,1) 
        plot(t, E); grid on; hold on; 
        plot(t, T); 
        plot(t, U); 
        ylabel('km^2/s^2')
        legend('Total', 'Kinetic', 'Potential') 
        title('Total Specific Energy: Kinetic - Potential') 
    subplot(2,1,2) 
        plot(t, [0 diff(E)]); grid on 
        title('Change in Total Specific Energy') 
        xlabel('Time (sec)') 
        ylabel('km^2/s^2')
sgtitle(name) 

save_pdf(h, 'prob5_energy')
    
%% subfunctions 

function save_pdf(h, name) 

% save as cropped pdf 
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,name,'-dpdf','-r0')
    
end 

