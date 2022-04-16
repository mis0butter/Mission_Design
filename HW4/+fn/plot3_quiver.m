function h = plot3_quiver(r1, r2, style, linew)

if ~exist('linew', 'var') 
    linew = 1; 
end 

if ~exist('style', 'var')
    h = quiver3(r1(1), r1(2), r1(3), r2(1), r2(2), r2(3), 'linewidth', linew);     
else 
    h = quiver3(r1(1), r1(2), r1(3), r2(1), r2(2), r2(3), style, 'linewidth', linew);     
end 

end 