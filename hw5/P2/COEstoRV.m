function outState = COEstoRV(A, MU)
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
%   A - 6x1 vector including:
%
%                   Semimajor Axis
%                   Eccentricity
%                   Inclination
%                   Right Ascention of the Ascending Node
%                   Argument of Periapse
%                   True Anomaly
%
%   MU - scalar gravitation parameter (G*M)
% 
% Outputs:
% 
%   R  - 3x1 vector representing position
%   V  - 3x1 vector representing velocity
%
% Assumptions/References:
%
%   NONE
%
% Dependencies:
%
%   rotation1()
%   rotation3()
%
% Modification History:
% 
%   17mar08     Brandon A. Jones      original version (header added)

   semi = A(1);
   e    = A(2);
   i    = A(3);
   node = A(4);
   arg  = A(5);
   true = A(6);
   
   p = semi*(1-e^2);  % p = semi-latus rectum
   
   RPQW(1) = p*cos(true) / (1+e*cos(true));
   RPQW(2) = p*sin(true) / (1+e*cos(true));
   RPQW(3) = 0;

   VPQW(1) = -sqrt(MU/p) * sin(true);
   VPQW(2) =  sqrt(MU/p) * (e+cos(true));
   VPQW(3) =  0;
   
   rot_matrix = rotation3(-node)*rotation1(-i)*rotation3(-arg);

   RIJK = rot_matrix*RPQW';
   VIJK = rot_matrix*VPQW';
   
   outState = [ RIJK(:); VIJK(:) ];
