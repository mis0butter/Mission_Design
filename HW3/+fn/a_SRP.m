function a_srp = a_SRP(et, X) 
% from Vallado 
% Special Perturbation Techniques Equation 8-45 

global const 

% get states --> Earth to Sun 
target   = 'Sun';
frame    = 'J2000';
observer = 'Earth';
abcorr   = 'NONE';
X_Esun   = fn.spice_state(et, target, frame, abcorr, observer); 
X_Esun   = X_Esun'; 

r_Esat = X(1:3); 
r_Esun = X_Esun(1:3); 

% Check if X is numeric or symbolic, then check if LOS is in Earth's shadow  
if isnumeric(X)
    % computing numerical - LOS depends 
    tau_min = (norm(r_Esat)^2 - dot(r_Esat, r_Esun)) / ( norm(r_Esat)^2 + norm(r_Esun)^2 - 2*dot(r_Esat, r_Esun) ); 
    LOS = false; 
    if tau_min < 0 || tau_min > 1
        LOS = true; 
    else
        ctau2 = ( (1-tau_min)*norm(r_Esat)^2 + dot(r_Esat, r_Esun)*tau_min ); 
        if ctau2 >= const.RE^2 
            LOS = true; 
        end 
    end 
else 
    % computing symbolic - LOS true 
    LOS = true; 
end 

r_Esat   = X(1:3); 
r_satsun = X_Esun(1:3) - r_Esat;

% normalize 
rhat_Esat   = r_Esat / norm(r_Esat); 
rhat_satsun = r_satsun / norm(r_satsun); 

% satellite coord unit vectors in ECI frame 
sat_z = rhat_Esat;                  % radial: minus-z is NADIR pointing, so plus-z is ZENITH pointing, or in r_Esat direction 
sat_vhat = X(4:6) / norm(X(4:6));   % direction of velocity is in orbital plane 
sat_y = cross(sat_z, sat_vhat);     % cross-track: perpendicular to orbital plane 
sat_x = cross(sat_y, sat_z);        % in-track: in orbital plane, in direction of motion 

% Get angle between sat x and z unit vectors and r_satsun vector 
r_satE  = -rhat_Esat; 
theta_z = acos(dot( sat_z, rhat_satsun ));  
theta_y = acos(dot( sat_y, rhat_satsun )); 
theta_x = acos(dot( sat_x, rhat_satsun )); 

% CHECK UNITS. USE KM !!!
theta_inc = 0; 
a_srp     = 0; 
a_srp_sp  = 0; 

% m_SC = 2000; 

if LOS 

    % solar panel SRP 
%     nhat = (X_Esun(1:3) - X(1:3))/norm(X_Esun(1:3) - X(1:3));  % normal to solar panel 
%     shat = X_Esun(1:3)/norm(X_Esun(1:3));                      % sun vector 
%     a_srp_sp = -const.p_srp * const.Asp/const.m_SC * cos(theta_inc) * (2*( const.Solar_Cells_Cd/3 + const.Solar_Cells_Cs*cos(theta_inc) )*nhat + (1-const.Solar_Cells_Cs)*shat); 
    a_srp_sp = calc_refl_srp( const.p_srp, const.Asp, theta_inc, const.m_SC, const.Solar_Cells_Cd, const.Solar_Cells_Cs, rhat_satsun, rhat_satsun ); 

end 
a_srp = a_srp + a_srp_sp; 

%% account for face & reflectivity 

if isnumeric(X)
    
    if LOS 
        
        % get SRPs 
        if theta_z < pi/2   % plus-z: White Paint  
            c_s  = const.White_Paint_Cs; 
            c_d  = const.White_Paint_Cd; 
            nhat = sat_z; 
            theta = theta_z; 
        else                % minus-z : Germanium Kapton 
            c_s  = const.Germanium_Kapton_Cs; 
            c_d  = const.Germanium_Kapton_Cd; 
            nhat = -sat_z; 
            theta = pi - theta_z; 
        end
        a_srp_z = calc_refl_srp( const.p_srp, const.Az, theta, const.m_SC, c_d, c_s, nhat, rhat_satsun ); 
        
        if theta_y < pi/2   % plus-y: MLI Kapton 
            nhat = sat_y; 
            theta = theta_y; 
        else                % minus-y: MLI Kapton 
            nhat = -sat_y; 
            theta = pi - theta_y; 
        end
        c_s = const.MLI_Kapton_Cs; 
        c_d = const.MLI_Kapton_Cd; 
        a_srp_y = calc_refl_srp( const.p_srp, const.Ay, theta, const.m_SC, c_d, c_s, nhat, rhat_satsun ); 
        
        if theta_x < pi/2   % plus-x: MLI Kapton
            nhat = sat_x; 
            theta = theta_x; 
        else                % minus-x: MLI Kapton 
            nhat = -sat_x; 
            theta = pi - theta_x; 
        end 
        c_s = const.MLI_Kapton_Cs; 
        c_d = const.MLI_Kapton_Cd; 
        a_srp_x = calc_refl_srp( const.p_srp, const.Ax, theta, const.m_SC, c_d, c_s, nhat, rhat_satsun ); 
        
        % Add SRP 
        a_srp = a_srp + a_srp_x + a_srp_y + a_srp_z; 
        
    end
    
end

a_srp = a_srp / 1000; 

end 

%% subfunctions 

function a_srp = calc_refl_srp( p_srp, A, theta, mass, c_Rd, c_Rs, nhat, shat )

a_srp = - p_srp * A * cos(theta) / mass * ... 
    ( 2 * ( c_Rd/3 + c_Rs * cos(theta) ) * nhat + (1 - c_Rs) * shat ); 

end

















