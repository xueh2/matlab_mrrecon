
function h = plotGTBullsEyeContoursOnMap(fmap, C, h, LineWidth, LineColor)
% h = plotGTBullsEyeContoursOnMap(fmap, C, LineWidth, LineColor)

if(nargin<3)
    h = 0;
end

if(nargin < 4)
    LineWidth = 1;
end

if(nargin < 5)
    LineColor = [1 0 0];
end

if(~ishandle(h))
    h = figure;
    imagescn(fmap, [], [], 12);
    hold on
end

for ii=1:numel(C)    
    
    starting_points = C(ii).starting_points;
    ending_points = C(ii).ending_points;
    
    l_ind = ii;
    if(l_ind>size(LineColor, 1))
        l_ind = size(LineColor, 1);
    end
    
    l_ind2 = ii;
    if(l_ind2>numel(LineWidth))
        l_ind2 = numel(LineWidth);
    end
    
    for tt=1:size(starting_points, 1)-1
       line([starting_points(tt,1) ending_points(tt,1)],[starting_points(tt,2) ending_points(tt,2)], ...
           'LineWidth', LineWidth(l_ind2), 'Color', LineColor(l_ind, :)); 
    end    
end

hold off
