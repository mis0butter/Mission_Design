function a = g_J2J3J4(X)

% geopotential due to J2 J3 J4 

    global const 
    
    r = sqrt( X(1)^2 + X(2)^2 + X(3)^2 ); 

    a0 = -const.muN;
    a2 = -3*const.J2*const.RN^2 / 2; 
    a3 = -const.J3*const.RN^3 / 2; 
    a4 = -5*const.J4*const.RN^4 / 8; 

    x1 = X(1); 
    x2 = X(2); 
    x3 = X(3); 

    dU1 = a0*x1/r^3 * ( 1 + a2/r^2*( 5*x3^2/r^2 - 1 ) + a3/r^3*( 35*x3^3/r^3 - 15*x3/r ) + a4/r^4*( 63*x3^4/r^4 - 42*x3^2/r^2 + 3 ) );  
    dU2 = a0*x2/r^3 * ( 1 + a2/r^2*( 5*x3^2/r^2 - 1 ) + a3/r^3*( 35*x3^3/r^3 - 15*x3/r ) + a4/r^4*( 63*x3^4/r^4 - 42*x3^2/r^2 + 3 ) ); 
    dU3 = a0*x3/r^3 * ( 1 + a2/r^2*( 5*x3^2/r^2 - 3 ) + a3/r^3*( 35*x3^3/r^3 - 30*x3/r + 3*r/x3 ) + a4/r^4*( 63*x3^4/r^4 - 70*x3^2/r^2 + 15 ));

    a = [dU1; dU2; dU3]; 
    
end 

%% 
