function [phi, h] = cylindrical2geodetic(rho, z, a, f, inDegrees)
%cylindrical2geodetic Geocentric cylindrical to geodetic coordinates
%
%       FOR INTERNAL USE ONLY -- This function is intentionally
%       undocumented and is intended for use only within other toolbox
%       functions and classes. Its behavior may change, or the function
%       itself may be removed in a future release.
%
%   [phi,h] = map.geodesy.internal.cylindrical2geodetic(rho,z,a,f,inDegrees)
%   returns geodetic coordinates corresponding to radial distance rho and
%   signed distance from the equator z in a spheroid-centric (ECEF)
%   cylindrical coordinate system.
%
%   Input Arguments
%   ----------------
%   rho -- Radial distance from polar axis of one or more points, specified
%      as a scalar value, vector, matrix, or N-D array. Values must be in
%      units that match the length unit of the semimajor axis.
%
%      Data types: single | double
%
%   z -- Signed distance from equatorial plane of one or more points from
%      the equatorial plane, specified as a scalar value, vector, matrix,
%      or N-D array (equivalent to the z-coordinate the spheroid-centric
%      ECEF Cartesian system). Values must be in units that match the
%      length unit of the semimajor axis.
%
%      Data types: single | double
%
%   a -- Semimajor axis of reference spheroid, specified as a scalar number.
%
%      Data type: double
%
%   f -- Flattening of reference spheroid, specified as a scalar number.
%
%      Data type: double
%
%   inDegrees -- Unit of angle flag, specified as a scalar logical. The
%      value true indicates that geodetic latitude phi is in degrees;
%      false indicates that phi is in radians.
%
%      Data type: logical
%
%   Output Arguments
%   ---------------
%   phi -- Geodetic latitude of one or more points, returned as a scalar
%      value, matrix, or N-D array. Units are determined by the inDegrees
%      flag. When in degrees, they lie in the closed interval [-90 90].
%
%   h -- Ellipsoidal height of one or more points, specified as a scalar
%      value, vector, matrix, or N-D array. Units are determined by the
%      length unit of the semimajor axis.
%
%   Notes
%   -----
%   This function follows standard elementwise behavior with respect to
%   inputs rho and z, including scalar expansion.
%
%   Longitude in a 3-D spheroid-centric cylindrical system is the same as
%   in the corresponding geodetic system, and hence is not needed as either
%   an input or an output. Another perspective is that this function
%   performs a 3-D spheroid-centric ECEF to geodetic transformation in the
%   plane of a meridian.
%
%   See also map.geodesy.internal.geodetic2cylindrical

% Copyright 2012-2019 The MathWorks, Inc.

%#codegen

if inDegrees
    sinfun = @sind;
    cosfun = @cosd;
    atan2fun = @atan2d;
else
    sinfun = @sin;
    cosfun = @cos;
    atan2fun = @atan2;
end

% Spheroid properties
b = (1 - f) * a;       % Semiminor axis
e2 = f * (2 - f);      % Square of (first) eccentricity
ep2 = e2 / (1 - e2);   % Square of second eccentricity

% Bowring's formula for initial parametric (beta) and geodetic
% (phi) latitudes
beta = atan2fun(z, (1 - f) * rho);
phi = atan2fun(z   + b * ep2 * sinfun(beta).^3,...
    rho - a * e2  * cosfun(beta).^3);

% Fixed-point iteration with Bowring's formula
% (typically converges within two or three iterations)
betaNew = atan2fun((1 - f)*sinfun(phi), cosfun(phi));
count = 0;
while any(beta(:) ~= betaNew(:)) && count < 5
    beta = betaNew;
    phi = atan2fun(z   + b * ep2 * sinfun(beta).^3,...
        rho - a * e2  * cosfun(beta).^3);
    betaNew = atan2fun((1 - f)*sinfun(phi), cosfun(phi));
    count = count + 1;
end

% Ellipsoidal height from final value for latitude
sinphi = sinfun(phi);
N = a ./ sqrt(1 - e2 * sinphi.^2);
h = rho .* cosfun(phi) + (z + e2 * N .* sinphi) .* sinphi - N;

% Implementation Notes from Rob Comer
% --------------------------------------------------
% The implementation below follows Wolf and DeWitt quite literally, with a
% few important exceptions required to ensure good numerical behavior:
%
% 1) I used ATAN2 (or ATAN2D) rather than ATAN in the formulas for beta and
%    phi. This avoids division by zero (or a very small number) for points
%    on (or near) the Z-axis.
%
% 2) Likewise, I used ATAN2 (or ATAN2D) instead of ATAN when computing beta
%    from phi (conversion from geodetic to parametric latitude), ensuring
%    stability even for points at very high latitudes.
%
% 3) Finally, I avoided dividing by cos(phi) -- also problematic at high
%    latitudes -- in the calculation of h, the height above the ellipsoid.
%    Wolf and Dewitt give
%
%                   h = sqrt(X^2 + Y^2)/cos(phi) - N,
%
%    or
%
%                   h = rho/cos(phi) - N,
%
%    The trick is to notice an alternative formula that involves
%    division by sin(phi) instead of cos(phi), then take a linear
%    combination of the two formulas weighted by cos(phi)^2 and
%    sin(phi)^2, respectively. This eliminates all divisions and,
%    because of the identity cos(phi)^2 + sin(phi)^2 = 1 and the
%    fact that both formulas give the same h, the linear
%    combination is also equal to h.
%
%    To obtain the alternative formula, we simply rearrange
%
%              Z = [N(1 - e^2) + h]sin(phi)
%    into
%              h = Z/sin(phi) - N(1 - e^2).
%
%    The linear combination is thus
%
%        h = (rho/cos(phi) - N) cos^2(phi)
%            + (Z/sin(phi) - N(1 - e^2))sin^2(phi)
%
%    which simplifies to
%
%        h = rho cos(phi) + Zsin(phi) - N(1 - e^2sin^2(phi)).
%
%    From here it's not hard to verify that along the Z-axis we
%    have h = Z - b and in the equatorial plane we have h = rho - a.