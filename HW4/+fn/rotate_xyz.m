function rf = rotate_xyz(ri, theta, axis)

if ~iscolumn(ri)
    ri = ri'; 
end 

if ~exist('axis', 'var')
    disp('Must give axis, no rotation')
    C = eye(3);     
end 

if axis == 1
    C = [   1   0           0; 
            0   cos(theta)  -sin(theta); 
            0   sin(theta)  cos(theta)  ]; 
elseif axis == 2
    C = [   cos(theta)  0   sin(theta); 
            0           1   0; 
            -sin(theta) 0   cos(theta)  ]; 
elseif axis == 3
    C = [   cos(theta) -sin(theta)  0; 
            sin(theta)  cos(theta)  0;
            0           0           1   ];
else
    disp('Must give axis = 1, 2, or 3. No rotation')
    C = eye(3); 
end 
    
rf = C * ri; 

end 