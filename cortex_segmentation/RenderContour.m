
function f = RenderContour(I, starting_points, ending_points, contourColor, contourWidth, f)
% overlay the contours on the anatomy slices
    
if ( exist('f') )
    figure(f);
else
    f=figure;
end
% f=gcf;
% f = 0;
minI = min(I(:));
maxI = 0.8*max(I(:));
if ( iscell(starting_points) == 0 )
    if (nargin < 4)
        contourColor = [1 0 0];
    end

    if (nargin < 5)
        contourWidth = 2;
    end

    num = size(starting_points, 1);
    hold on
    if (isempty(I) == 0)
        imshow(I,[minI maxI]);
    end
    for tt=1:num
       line([starting_points(tt,1) ending_points(tt,1)],[starting_points(tt,2) ending_points(tt,2)], ...
           'LineWidth', contourWidth, 'Color', contourColor); 
    end
    
%     f = line([starting_points(:,1) ending_points(:,1)],[starting_points(:,2) ending_points(:,2)], ...
%            'LineWidth', contourWidth, 'Color', contourColor); 
    
    hold off
    axis image
%     drawnow;
    return
else
    
    num = numel(starting_points);
    
    hold on
    if (isempty(I) == 0)
        imshow(I,[minI maxI]);
    end

    for ss=1:num 
        ps = starting_points{ss};
        pe = ending_points{ss};
        
        N = size(ps, 1);
        for tt=1:N
           line([ps(tt,1) pe(tt,1)],[ps(tt,2) pe(tt,2)], ...
               'LineWidth', contourWidth(ss), 'Color', contourColor(ss,:)); 
        end
    end
%     drawnow;
    hold off
end
axis image
return;