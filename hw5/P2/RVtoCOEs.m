function [COEs] = RVtoCOEs(R,V,MU)
%
% function [COEs] = RVtoCOEs(R,V,MU)
% ---------------------------------------------------------------------
% 
% Description:
%
%   Function to convert the provide radius and velocity vector to classical
%   orbital elements (COEs).
% 
% Inputs:
%
%   R  - 3x1 vector representing position
%   V  - 3x1 vector representing velocity
%   MU - scalar gravitation parameter (G*M)
% 
% Outputs:
% 
%   A - 6x1 vector including:
%
%                   Semimajor Axis
%                   Eccentricity
%                   Inclination
%                   Right Ascention of the Ascending Node
%                   Argument of Periapse
%                   True Anomaly
%
% Assumptions/References:
%
%   NONE
%
% Dependencies:
%
%   NONE
%
% Modification History:
% 
%   17mar08     Brandon A. Jones      original version (header added)


%  Just make sure we have column vectors for future processing.
R = R(:);
V = V(:);
dotRV = R'*V;

%  Precompute the radius and speed
radius = norm(R);
speed  = norm(V);

%  Get the energy and the SMA
energy=speed^2/2 - MU/radius;
a = -MU/2/energy;

%  Get the angular momentum vector
H = cross(R,V);

%  Get the eccentrcity vectory
E = ((speed*speed-MU/radius)*R - (dotRV)*V)/MU';
ecc = norm(E);

%  Get the inclination
inclination = arccos(H(3)/norm(H));

%  Now, we get the RAAN
k = [0,0,1];
N = cross(k,H);
normN = norm(N);
node = arccos(N(1)/normN);
if N(2)<0
    node=2*pi-node;
end;


%  Get the argument of perigee
arg = arccos(  N(:)'*E(:) / ( normN*ecc )  );
if E(3)<0
    arg = 2*pi - arg;
end;

%  Get the true anomaly
true = arccos( (E(:)'*R(:))/(ecc*radius) );
if dotRV<0
    true = 2*pi - true;
end

%  Now, we need to check for the singular cases
if ecc < 1e-10 && inclination < 1e-10
    
    %  Set these values to zero.
    node = 0.0;
    arg  = 0.0;
    
    %  Compute the true longitude of periapsis
    true = arccos( R(1)/radius );
    if R(2) < 0
        true = 2*pi - true;
    end
    
elseif ecc < 1e-10
    
    %  Set the singular elements to zero
    arg = 0.0;
    
    %  Compute the argument of latitude
    true = arccos( dot(N,R)/(normN*radius) );
    if R(3) < 0
        true = 2*pi - true;
    end
    
elseif inclination < 1e-10
    
    %  Set the singular elements to zero
    node = 0.0;
    
    %  Compute the Longitude of periapsis
    arg = arccos( E(1)/ecc );
    if E(2) < 0
        arg = 2*pi - arg;
    end
    
end

COEs = [a; ...
        ecc; ...
        inclination; ...
        node; ...
        arg; ...
        true ];
    
end
    
    function retVal = arccos( argument )
        if abs(argument-1) < 1e-12
            argument = 1;
        elseif abs(argument+1) < 1e-12
            argument = -1;
        end
        retVal = acos(argument);
    end
