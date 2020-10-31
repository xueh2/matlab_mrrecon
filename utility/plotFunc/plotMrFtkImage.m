
function [xyAxes, xyImageHandle, endo_handle, epi_handle, landmark_handle] = plotMrFtkImage(im, header, h, showAxis, endo_pt, epi_pt, landmark_pt)

if(nargin<5)
    endo_pt = [];
end

if(nargin<6)
    epi_pt = [];
end

if(nargin<7)
    landmark_pt = [];
end

endo_handle = -1;
epi_handle = -1;
landmark_handle = -1;

xlim3D(1) = header.spacingX*-0.5;
xlim3D(2) = header.spacingX*(header.sizeX-0.5);

ylim3D(1) = header.spacingY*-0.5;
ylim3D(2) = header.spacingY*(header.sizeY-0.5);

minLim = min([xlim3D(1) ylim3D(1)]);
maxLim = max([xlim3D(2) ylim3D(2)]);

maxXYLim = max([ylim3D(2)-ylim3D(1) xlim3D(2)-xlim3D(1)]);

if ( ~ishandle(h) )
    figure;
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0.02 0.05 0.95 0.95]);
    xyAxes = gca;
else
    xyAxes = h;
end
% xy
axes(xyAxes); colormap(gray);
hold on
axis(xyAxes, 'ij', 'equal');
if ( ~showAxis )
    axis(xyAxes, 'off');
    set(xyAxes, 'Position', [0 0 1.0 1.0]);
else
    set(xyAxes, 'Position', [0.1 0.1 0.8 0.8]);
end

box(xyAxes, 'off');
set(xyAxes, 'View', [0 90]);
set(xyAxes, 'XAxisLocation', 'top', 'YAxisLocation', 'left'); 

% set(xyAxes, 'XLim', [minLim maxLim] );
% set(xyAxes, 'YLim', [minLim maxLim] );
set(xyAxes, 'XLim', xlim3D );
set(xyAxes, 'YLim', ylim3D );
set(xyAxes, 'XColor', [1 0 0], 'YColor', [0 0 1]);

vi = im;
x = [xlim3D(1) xlim3D(2); xlim3D(1) xlim3D(2)];
y = [ylim3D(1) ylim3D(1); ylim3D(2) ylim3D(2)];
axes(xyAxes);
xyImageHandle = surface('XData',x,'YData',y,'ZData',[0 0; 0 0],'CData', double(vi),'FaceColor','texturemap','EdgeColor','none');

if(~isempty(endo_pt))
    endo_handle = plot(xyAxes, endo_pt(:,1), endo_pt(:,2), 'r');
end

if(~isempty(epi_pt))
    epi_handle = plot(xyAxes, epi_pt(:,1), epi_pt(:,2), 'g');
end

if(~isempty(landmark_pt))
    landmark_handle = plot(xyAxes, landmark_pt(:,1), landmark_pt(:,2), 'r+', 'MarkerSize', 14);
end

hold off

set(gcf, 'Renderer', 'OpenGL');
