function h = plotBoundingBoxOf3DVolume(data, header, h)
%============================
% h = plotBoundingBoxOf3DVolume(data, header, h)
% This is a function for drawing bounding box a 3D volume
% information.
% Input:
%   data : [column row SLC], [2nd dimension, 1st dimension, 3rd dimension]
%   header : the ftk header format, including positionPatient (tx ty tz) and orientationPatient (row vector; col vector; slc vector)
%   figh = figure handle
%   Note: all access starts from 0
%======================================================================

if ( ~ishandle(h) )
    figure;
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0.02 0.05 0.95 0.95], 'Renderer', 'zbuffer');
    h = gca;
    
    axis(h, 'ij', 'equal');
    set(h, 'View', [45 45]);
    set(h, 'XAxisLocation', 'top', 'YAxisLocation', 'left'); 

    set(h, 'XColor', [1 0 0], 'YColor', [0 0 1], 'ZColor', [0 0 0]);    
else
    set(gcf, 'Renderer', 'zbuffer');
end

[wc1, wc2, wc3, wc4, wc5, wc6, wc7, wc8] = computeEightCorner(data, header);

hold on

x = [wc1(1) wc2(1) wc3(1) wc4(1) wc1(1)];
y = [wc1(2) wc2(2) wc3(2) wc4(2) wc1(2)];
z = [wc1(3) wc2(3) wc3(3) wc4(3) wc1(3)];
line(x, y, z, 'LineStyle', '-', 'LineWidth', 2.0);

x = [wc5(1) wc6(1) wc7(1) wc8(1) wc5(1)];
y = [wc5(2) wc6(2) wc7(2) wc8(2) wc5(2)];
z = [wc5(3) wc6(3) wc7(3) wc8(3) wc5(3)];
line(x, y, z, 'LineStyle', '-', 'LineWidth', 2.0);

x = [wc2(1) wc6(1) wc7(1) wc3(1) wc2(1)];
y = [wc2(2) wc6(2) wc7(2) wc3(2) wc2(2)];
z = [wc2(3) wc6(3) wc7(3) wc3(3) wc2(3)];
line(x, y, z, 'LineStyle', '-', 'LineWidth', 2.0);

x = [wc1(1) wc4(1) wc8(1) wc5(1) wc1(1)];
y = [wc1(2) wc4(2) wc8(2) wc5(2) wc1(2)];
z = [wc1(3) wc4(3) wc8(3) wc5(3) wc1(3)];
line(x, y, z, 'LineStyle', '-', 'LineWidth', 2.0);

hold off

return

