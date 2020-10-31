
function f = RenderContour2(I, starting_points, ending_points, contourColor, contourWidth, pts_Image, pt_color, pt_size, pt_symbol, slicenum, f)
% overlay the contours on the anatomy slices
    
if ( exist('f') & ~isempty(f) & ishandle(f) )
    figure(f);
else
    f=figure;
end
[H, W] = size(I);

% set(f, 'MenuBar', 'none' );

set(f, 'Units', 'pixels');
set(f, 'position', [200  200  1.1*W  1.1*H]);

hold on   
h = get(f, 'CurrentAxes');
set(h, 'Units', 'pixels');
set(h, 'position',[5  5  W  H]);

if (isempty(I) == 0)
    imshow(uint8(I));
end
    
if ( iscell(starting_points) == 0 )
    if (nargin < 4)
        contourColor = [1 0 0];
    end

    if (nargin < 5)
        contourWidth = 2;
    end

    num = size(starting_points, 1);
    
%     if (isempty(I) == 0)
%         imshow(uint8(I));
%     end
    for tt=1:num
       line([starting_points(tt,1) ending_points(tt,1)],[starting_points(tt,2) ending_points(tt,2)], ...
           'LineWidth', contourWidth, 'Color', contourColor); 
    end
    
else
    
    num = numel(starting_points);
    
%     hold on
%     h = get(f, 'CurrentAxes');
%     set(h, 'Units', 'pixels');
%     set(h, 'position',[5  5  W  H]);

    for ss=1:num 
        ps = starting_points{ss};
        pe = ending_points{ss};
        
        N = size(ps, 1);
        for tt=1:N
           line([ps(tt,1) pe(tt,1)],[ps(tt,2) pe(tt,2)], ...
               'LineWidth', contourWidth(ss), 'Color', contourColor(ss,:)); 
        end
    end
    
%     x = pts_Image(:, 1);
%     y = pts_Image(:, 2);
%     plot(x, y,  'MarkerEdgeColor', pt_color, 'MarkerFaceColor', pt_color, ...
%                 'MarkerSize', pt_size, 'LineStyle', 'none', 'Marker', pt_symbol)
end

if ( isempty(pts_Image)==0 )
    x = pts_Image(:, 1);
    y = pts_Image(:, 2);
    plot(x, y,  'MarkerEdgeColor', pt_color, 'MarkerFaceColor', pt_color, ...
                    'MarkerSize', pt_size, 'LineStyle', 'none', 'Marker', pt_symbol)
end    
if ( slicenum > 0 )
    
    text(0.1*W, 0.9*H, num2str(slicenum), 'FontSize', 14, 'Color', [1 1 1]);
    
end
hold off
axis image
axis ij
return;