function y = G_bias_fn(X, XS, STA)
% X  = satellite coords (SAME COORD)
% XS = station coords (SAME COORD)
 
KJL_rbias = 0;  
DGO_rbias = 0; 
% ACB_rbias = 0; 
ACB_rbias = 0.020; 

if      STA == 1;   rbias = KJL_rbias;   
elseif  STA == 2;   rbias = DGO_rbias;   
elseif  STA == 3;   rbias = ACB_rbias;   
end

r_site = [X(1)-XS(1); X(2)-XS(2); X(3)-XS(3)]; 
v_site = [X(4)-XS(4); X(5)-XS(5); X(6)-XS(6)]; 
d      = norm(r_site) - rbias; 
v      = dot(v_site, r_site / d); 

y      = [d; v]; 

end 
    
    