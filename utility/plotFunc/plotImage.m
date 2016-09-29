
function plotMrFtkImage(im, header)

xlim3D(1) = header.spacingX*-0.5;
xlim3D(2) = header.spacingX*(header.sizeX-0.5);

ylim3D(1) = header.spacingY*-0.5;
ylim3D(2) = header.spacingY*(header.sizeY-0.5);

minLim = min([xlim3D(1) ylim3D(1)]);
maxLim = max([xlim3D(2) ylim3D(2) zlim3D(2)]);

maxXYLim = max([ylim3D(2)-ylim3D(1) xlim3D(2)-xlim3D(1)]);

figure;
xyAxes = gca;
% xy
axes(xyAxes); colormap(gray);
axis(xyAxes, 'ij', 'equal');
box(xyAxes, 'on');
set(xyAxes, 'View', [0 90]);
set(xyAxes, 'XAxisLocation', 'top', 'YAxisLocation', 'left'); 

set(xyAxes, 'XLim', [minLim maxLim] );
set(xyAxes, 'YLim', [minLim maxLim] );
set(xyAxes, 'XColor', [1 0 0], 'YColor', [0 0 1]);

vi = im;
x = [xlim3D(1) xlim3D(2); xlim3D(1) xlim3D(2)];
y = [ylim3D(1) ylim3D(1); ylim3D(2) ylim3D(2)];
axes(xyAxes);
xyImageHandle = surface('XData',x,'YData',y,'ZData',[0 0; 0 0],'CData', double(vi),'FaceColor','texturemap','EdgeColor','none');

set(gcf, 'Renderer', 'OpenGL');
