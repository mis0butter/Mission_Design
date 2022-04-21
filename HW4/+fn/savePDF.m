function savePDF(h, plot_name, plot_path)
% ------------------------------------------------------------------------
% Crop and save MatLAB figure as PDF
% 
% Ex: 
%   >> fn.savePDF(gcf) 
% ------------------------------------------------------------------------

if ~exist('plot_name', 'var')
    plot_name = h.Name; 
end 

% check if directory exists, if not create one
if ~exist('plot_path', 'var')
    plot_path = ''; 
elseif ~exist('plot_path', 'dir')
    mkdir(plot_path); 
end

% figure handle 
fig = h; 
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'-dpdf','-painters','-r600','-bestfit',strcat(plot_path,plot_name,'.pdf'));
% print(fig,'-dpdf','-painters','-r600','-bestfit',strcat(plot_name));
end

% %% Crop and save MatLAB figure as PDF
% function savePDF(plot_path,plot_name)
% % check if directory exists, if not create one
% if ~exist(plot_path, 'dir')
%     mkdir(plot_path)
% end
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'-dpdf','-painters','-r600','-bestfit',strcat(plot_path,plot_name));
% end