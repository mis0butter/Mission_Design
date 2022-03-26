function [ell_1, ell_2] = bett_lambert(r_sunE, r_sunM, mu)

%% bettadpur lambert 

% Lambert - from VALLADO ed. 4, pgs 467 - 475 

% Departure (Earth)
rE_mag = norm(r_sunE); 

% Arrival (Mars) 
rM_mag = norm(r_sunM); 

% Vallado method ... 
cos_dv = dot(r_sunE, r_sunM) / (rE_mag * rM_mag); 

% chord 
c = sqrt( rE_mag^2 + rM_mag^2 - 2*rE_mag*rM_mag*cos_dv ); 

% semiperimeter 
s = ( rM_mag + rE_mag + c ) / 2; 

% min semimajor axis 
amin = s/2;  

% mean motion 
n  = sqrt(mu/amin^3); 

% time of flight 
alpha1 = 2 * asin( sqrt( s/(2*amin) ) ); 
alpha2 = 2*pi - alpha1; 
beta1  = 2 * asin( sqrt( (s-c)/(2*amin) ) ); 
beta2  = - beta1; 

dE1 = alpha1 - beta1; 
dE2 = alpha2 - beta2; 

% solution 4 
t_a1 = alpha1 - sin(alpha1); 
t_a2 = alpha2 - sin(alpha2); 
t_b1 = beta1 - sin(beta1); 
t_b2 = beta2 - sin(beta2); 

d_t1 = 1/n * (t_a1 + t_b1); 
d_t2 = 1/n * (t_a2 + t_b2); 

f1 = 1 - amin / rE_mag * ( 1 - cos(dE1) ); 
f2 = 1 - amin / rE_mag * ( 1 - cos(dE2) ); 
g1 = d_t1 - 1/n * ( dE1 - sin(dE1) );  
g2 = d_t2 - 1/n * ( dE2 - sin(dE2) );  

v_E1 = 1/g1 * ( r_sunM - f1 * r_sunE );
v_E2 = 1/g2 * ( r_sunM - f2 * r_sunE );

% ellipse 1 - rv is from Sun to satellite 
rv0_e1 = [r_sunE, v_E1]; 
rv0_e2 = [r_sunE, v_E2]; 

ell_1.rv0 = rv0_e1; 
ell_1.tof = d_t1; 
ell_2.rv0 = rv0_e2; 
ell_2.tof = d_t2; 


oe_e1 = rvOrb.rv2orb(rv0_e1', mu); 
oe_e2 = rvOrb.rv2orb(rv0_e2', mu); 

end 