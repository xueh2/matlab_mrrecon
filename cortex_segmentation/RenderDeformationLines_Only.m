
function f = RenderDeformationLines_Only(I, starting_points, ending_points, LineWidth, LineColor, pt_size, ...
    starting_pt_color,  starting_pt_symbol, ending_pt_color,  ending_pt_symbol, ...
    slicenum, f)
% overlay the contours on the anatomy slices
    
if ( exist('f') & ~isempty(f) & ishandle(f) )
    figure(f);
else
    f=figure;
end
[H, W] = size(I);

% set(f, 'MenuBar', 'none' );

% set(f, 'Units', 'pixels');
% set(f, 'position', [200  200  1.1*W  1.1*H]);

hold on   
% h = get(f, 'CurrentAxes');
% set(h, 'Units', 'pixels');
% set(h, 'position',[5  5  W  H]);

% if (isempty(I) == 0)
%     imshow(uint8(I));
% end
    
num = size(starting_points, 1);
    
for tt=1:num
        line([starting_points(tt,1) ending_points(tt,1)],[starting_points(tt,2) ending_points(tt,2)], ...
           'LineWidth', LineWidth, 'Color', LineColor);       
end
    
x = starting_points(:, 1);
y = starting_points(:, 2);
plot(x, y,  'MarkerEdgeColor', starting_pt_color, 'MarkerFaceColor', starting_pt_color, ...
                'MarkerSize', pt_size, 'LineStyle', 'none', 'Marker', starting_pt_symbol);

% x = ending_points(:, 1);
% y = ending_points(:, 2);
% plot(x, y,  'MarkerEdgeColor', ending_pt_color, 'MarkerFaceColor', ending_pt_color, ...
%                 'MarkerSize', pt_size, 'LineStyle', 'none', 'Marker', ending_pt_symbol);
    
% if ( slicenum > 0 )
%     
%     text(0.1*W, 0.9*H, num2str(slicenum), 'FontSize', 14, 'Color', [1 1 1]);
%     
% end
hold off
axis image

return;