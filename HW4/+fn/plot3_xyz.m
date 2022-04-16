function h = plot3_xyz(X, style, linew, i)
% ------------------------------------------------------------------------
% PURPOSE: 
%   Create 3D plot of position vector or time history of position vectors
%   ([N x 3] matrix) 
% 
% INPUTS: 
%   X       = [N x 3] position vectors 
%   style   = line style, e.g. 'b' or 'r--' 
%   linew   = linewidth, e.g. 1 
%   i       = index. If blank function plots all (:). if i is string 'end',
%               plots last index 
% 
% OUTPUTS: 
%   h       = line output handle 
% ------------------------------------------------------------------------

if ~exist('style', 'var') 
    style = ''; 
end 

if ~exist('linew', 'var')
    linew = 1; 
end 

if ~exist('i', 'var')
    h = plot3(X(:,1), X(:,2), X(:,3), style, 'linewidth', linew); 
elseif strcmp(i, 'end')
    h = plot3(X(end,1), X(end,2), X(end,3), style, 'linewidth', linew);     
else
    h = plot3(X(i,1), X(i,2), X(i,3), style, 'linewidth', linew); 
end 

end 