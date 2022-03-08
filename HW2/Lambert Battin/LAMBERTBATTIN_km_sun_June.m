function [vd, va] = LAMBERTBATTIN_km_sun_June(rd, ra, dm, dt)

% mu = 3.986004418e14;   % m3/s2
mu_sun_m = 1.32712440018e20; 
mu_sun_km = mu_sun_m / (1000^3); 
mu = mu_sun_km; 

% cos and sin dnu 
ra_mag  = norm(ra);
rd_mag  = norm(rd);
cos_dnu = dot(rd,ra)/(rd_mag*ra_mag);
dnu     = acos(cos_dnu); 
sin_dnu = sin(dnu); 

% the angle needs to be positive to work for the long way
if dnu < 0.0
    dnu = 2*pi + dnu;
end

% geometry 
c   = sqrt( rd_mag^2 + ra_mag^2 - 2*rd_mag*ra_mag*cos_dnu ); 
s   = ( rd_mag + ra_mag + c ) / 2; 
eps = (ra_mag - rd_mag) / rd_mag; 

% tan2w --> rdp 
sqrt_rda = sqrt( ra_mag / rd_mag ); 
tan22w    = ( eps^2/4 ) / ( sqrt_rda + sqrt_rda^2*( 2 + sqrt_rda ) ); 
rdp      = sqrt( rd_mag*ra_mag ) * ( cos( dnu/4 )^2 + tan22w ); 

% obtain L 
if dnu < pi 
    num = sin(dnu/4)^2 + tan22w; 
    den = sin(dnu/4)^2 + tan22w + cos(dnu/2); 
else
    num = cos(dnu/4)^2 + tan22w - cos(dnu/2); 
    den = cos(dnu/4)^2 + tan22w; 
end 
L = num / den; 

m = mu * dt^2 / ( 8*rdp^3 ); 

% orbit is elliptical 
xn = L;
eta = xn / ( sqrt( 1+xn ) + 1 )^2;  
x = -xn; 

i = 0; 
% loop 
while (1) 
    
    i = i + 1; 
    
    x = xn; 
    
    % omg crazy recursion in fraction denominator 
    tempx = seebatt(x); 

    % h1 
    num = (L+x)^2 * (1 + 3*x + tempx); 
    den = ( 1 + 2*x + L ) * ( 4*x + tempx*( 3+x ) ); 
    h1  = num / den; 
    
    % h2 
    num = m*(x - L + tempx); 
    h2  = num / den; 
    
    % solve cubic y^3 - y^2 - h1*y^2 - h2 = 0 
    rts = roots([1 -1 -h1 -h2]); 
    
    B = 27*h2 / ( 4*( 1 + h1 )^3 ); 
    U = B / ( 2*( sqrt( 1+B ) + 1 ) ); 
    
    K = seebattk(U); 
    y = (1 + h1)/3 * ( 2 + ( sqrt(1+B) )/( 1 + 2*U*K^2 ) ); 
    xn = sqrt( ( (1-L)/2 )^2 + ( m/y^2 ) ) - (1+L)/2; 
    
    if abs(xn - x) < 0.000001 && i > 30
        break 
    end 
        
end 

% semi-major axis!!! 
a = mu_sun_km*dt^2 / ( 16 * rdp^2 * x * y^2 ); 

if a > 0 
    
    % obtain f, g, and dg 
    alpha = 2 * asin( sqrt( s/(2*a) ) ); 
    beta  = 2 * asin( sqrt( (s-c)/(2*a) ) ); 
    
    % min 
    amin = s / 2; 
    tmin = sqrt( amin^3/mu_sun_km ) * ( pi - beta + sin(beta) ); 
    if dt > tmin 
        alpha = 2*pi - alpha; 
    end 
    de    = alpha - beta; 
    
    f  = 1 - a/rd_mag * ( 1 - cos(de) ); 
    g  = dt - sqrt( a^3/mu_sun_km ) * ( de - sin(de) ); 
    dg = 1 - a/ra_mag * ( 1 - cos(de) ); 

else
    
    disp('oh no') 
    
end 

% velocities!! 
vd = (ra - f*rd) / g; 
va = (dg*ra - rd) / g; 

end 