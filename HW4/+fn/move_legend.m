function move_legend 
% ------------------------------------------------------------------------ 
% Purpose: Move legend outside of plot area and resize all axes (subplots)
% to be the same. 
% ------------------------------------------------------------------------ 

    % Get all axes handles 
    h_axes = findobj(gcf, 'type', 'axes'); 
    
    for i = 1:length(h_axes) 
        
        % If the axes is NOT 'suptitle' 
        if ~isequal(get(h_axes(i), 'tag'), 'suptitle')
            
            % Get current axes handle 
            h = h_axes(i); 
            
            % Shrink subplot 
            h_pos       = get(h, 'position'); 
            new_h_pos   = [h_pos(1), h_pos(2), 0.85*h_pos(3), h_pos(4)]; 
            set(h, 'position', new_h_pos); 
            
            % If subplot contains a legend 
            if ~isempty(h.Legend)
                
                % Get and set new h_legend position 
                h_leg_pos        = get(h.Legend, 'position'); 
                new_leg_pos(1)   = new_h_pos(1) + new_h_pos(3) + 0.02*h_pos(3); 
                new_leg_pos(2:4) = h_leg_pos(2:4); 
                set(h.Legend, 'position', new_leg_pos); 
                
            end 
        end 
    end 
    
end 