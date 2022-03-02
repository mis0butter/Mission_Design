%% 2: Download the SPICE toolkit in your preferred language from https://naif.jpl.nasa.gov/naif/toolkit.html.
% (MATLAB version is ‘MICE’; what follows assume you will use that, but feel free to use any of the available languages/platforms).
% Also download these three data files from Canvas\files\spice\ (or from the web):
% 'de421.bsp'
% 'pck00010.tpc'
% 'naif0011.tls'
% Familiarize* yourself with (at least these) functions:
% cspice_furnsh()
% mice_spkezr()
% cspice_pxform()
% cspice_str2et()
% *A list of all functions is here in the documentation: mice/doc/html/mice/index.html
% Use 'NONE' for the aberration correction input. For ii) – iv), consider t0 = 'Oct 20, 2020 11:00 AM CST' and use axis
% equal and view(-40,30).

% CSPICE_FURNSH loads SPICE kernel files from FILE into MATLAB.

% MICE_SPKEZR returns the state (position and velocity) of
%    a target body relative to an observing body, optionally
%    corrected for light  time (planetary aberration) and stellar
%    aberration.

% CSPICE_STR2ET converts a string representing an epoch to a
%    double precision value representing the number of TDB seconds
%    past the J2000 epoch corresponding to the input epoch.

clear; clc;

%  Load kernel file 
cspice_furnsh( 'spice_data/naif0011.tls' )
cspice_furnsh( 'spice_data/de421.bsp' )       
cspice_furnsh( 'spice_data/pck00010.tpc ') 

%  Define parameters for a state lookup:
t0          = 'Oct 20, 2020 11:00 AM CST'; 
abcorr      = 'NONE';


%% 2i: Convert t0 to ET, i.e. seconds past J2000, the base time variable 
% for SPICE. function calls.

%  Convert the epoch to ephemeris time. 
et_t0       = cspice_str2et( t0 );


%% 2ii: Plot the trajectories of Mercury, Venus and Earth with respect to 
% the sun in the inertial 'ECLIPJ2000' frame from t=t0 to t0+700 days.

% Mercury 
target      = 'Mercury';
frame       = 'ECLIPJ2000';
observer    = 'Sun';
days        = 700; 

[matMS, stargMS, rvMS] = spice_state(days, et_t0, target, frame, abcorr, observer); 
rMSx = rvMS(:,1); rMSy = rvMS(:,2); rMSz = rvMS(:,3); 

% Venus 
target      = 'Venus';

[matVS, stargVS, rvVS] = spice_state(days, et_t0, target, frame, abcorr, observer); 
rVSx = rvVS(:,1); rVSy = rvVS(:,2); rVSz = rvVS(:,3); 

% Earth 
target      = 'Earth';

[matES, stargES, rvES] = spice_state(days, et_t0, target, frame, abcorr, observer); 
rESx = rvES(:,1); rESy = rvES(:,2); rESz = rvES(:,3); 

% ------------------------------------------------------------------------ 
% Plot 

fname = sprintf('Problem 2: %s trajectory wrt %s in %s frame', target, observer, frame); 
figure('name', fname); 
    plot3(rMSx, rMSy, rMSz, '-o', 'markersize', 2); 
    hold on; grid on; 
    plot3(rVSx, rVSy, rVSz, '-o', 'markersize', 2); 
    plot3(rESx, rESy, rESz, '-o', 'markersize', 2); 
    legend('Mercury', 'Venus', 'Earth', 'location', 'southoutside')
    
    view(-40, 30); axis equal; 
    title(fname) 
    xlabel('x (LU)'); ylabel('y (LU)'); zlabel('z (LU)'); 
    
    
%% 2iii: Plot the Earth’s position relative to the moon’s position in the 
% inertial 'J2000' frame from t=t0 to t0+27 days. Discuss results.

% Earth wrt Moon 
target      = 'Earth';
frame       = 'J2000';
observer    = 'Moon';
days        = 27; 

[matEM, stargEM, rvEM] = spice_state(days, et_t0, target, frame, abcorr, observer); 
rEMx = rvEM(:,1); rEMy = rvEM(:,2); rEMz = rvEM(:,3); 

% ------------------------------------------------------------------------ 
% Plot 

fname = sprintf('Problem 2: %s trajectory wrt %s in %s frame', target, observer, frame); 
figure('name', fname); 
    hold on; grid on; 
    plot3(rEMx, rEMy, rEMz, '-o', 'markersize', 4, 'linewidth', 1.5); 
    plot3(0, 0, 0, 'o', 'markersize', 4, 'linewidth', 1.5)
    legend('Earth', 'Moon', 'location', 'southoutside')
    
    view(-40, 30); axis equal; 
    title(fname) 
    xlabel('x (LU)'); ylabel('y (LU)'); zlabel('z (LU)'); 
    
    
%% 2iv: Plot the Earth’s position relative to the moon’s position in the 
% body fixed 'IAU_Moon' frame from t=t0 to t0+27 days. use axis equal; 
% view(-40,30). Discuss results.

% Earth wrt Moon 
target      = 'Earth';
frame       = 'IAU_Moon';
observer    = 'Moon';
days        = 27; 

[matEM, stargEM, rvEM] = spice_state(days, et_t0, target, frame, abcorr, observer); 
rEMx = rvEM(:,1); rEMy = rvEM(:,2); rEMz = rvEM(:,3); 

% ------------------------------------------------------------------------ 
% Plot 

fname = sprintf('Problem 2: %s trajectory wrt %s in IAU Moon frame', target, observer); 
figure('name', fname); 
    hold on; grid on; 
    plot3(rEMx, rEMy, rEMz, '-o', 'markersize', 4, 'linewidth', 2); 
%     plot3(0, 0, 0, 'o', 'markersize', 4, 'linewidth', 1.5)
%     legend('Earth', 'Moon', 'location', 'southoutside')
    
    view(-40, 30); axis equal; 
    title(fname) 
    xlabel('x (LU)'); ylabel('y (LU)'); zlabel('z (LU)'); 

    
%% 2v: Use the code to verify the apparent retrograde path of Mars across 
% the night sky: https://mars.nasa.gov/all-aboutmars/night-sky/retrograde/. 
% Approximately reproduce the 2005 chart illustrated here: ... 

% Mars wrt Earth 
target      = 'Mars';
frame       = 'J2000';
observer    = 'Earth';
days        = 200; 

% Convert the epoch to ephemeris time.
t0      = 'July 27, 2005'; 
et_t0   = cspice_str2et( t0 );
% July 27, 2005 + 200 days = Feb 12, 2006 

[matME, stargME, rvME] = spice_state(days, et_t0, target, frame, abcorr, observer); 
rMEx = rvME(:,1); rMEy = rvME(:,2); rMEz = rvME(:,3); 

for i = 1:length(rvME) 
    
    RA(i,:) = atand(rMEy(i,:)/rMEx(i,:)); 

%     If Xq is negative then add 180 degrees to alpha  
%     If Xq is positive and Yq is negative then add 360 degrees to
%     alpha
% 
%     alpha is usually expressed in hours, so divide by 15

    dec(i,:) = atand( rMEz(i,:) / sqrt(rMEx(i,:)^2 + rMEy(i,:)^2)); 

end 

% ------------------------------------------------------------------------ 
% Plot 

fname = sprintf('Problem 2: %s retrograde in %s night sky', target, observer); 
figure('name', fname); 
    plot(RA, dec); 
    hold on; grid on; 
    plot(RA(1), dec(1), 'o', 'linewidth', 3); 
    plot(RA(end), dec(end), 'x', 'markersize', 10, 'linewidth', 3); 
    
    title(fname) 
    legend('Mars retrograde path', 'July 27, 2005', 'Feb 12, 2006', ... 
        'location', 'northwest'); 
    xlabel('RA (deg)')
    ylabel('Declination (deg)'); 
    

%% functions 

function [mat, starg, rv] = spice_state ... 
    (duration, et_t0, target, frame, abcorr, observer)

    for i = 1:duration + 1

        % Propagate 1 day forward in time 
        epoch       = et_t0 + (i-1) * 24*60*60; 

        % CSPICE_PXFORM returns the matrix that transforms position
        % vectors from one specified frame to another at a specified epoch.
        % Retrieve the transformation matrix FROM (1) TO (2) at epoch (3).
        mat{i,:}    = cspice_pxform( 'IAU_MERCURY', 'ECLIPJ2000', epoch );

        %  Look-up the state for the defined parameters.
        starg{i,:}  = mice_spkezr( target, epoch, frame, abcorr, observer);
        rv(i,:)      = starg{i,:}.state(1:6); 

    end 

end
