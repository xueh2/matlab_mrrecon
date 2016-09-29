function h = plotMrFtkImageIn3D(data, header, h, centre, width)
%============================
% [ fighandle, surfh ] = DrawImageIn3D(data, header, fighandle, transparency)
% This is a function for drawing slices in 3D space based on the dicom coordinates
% information.
% Input:
%   data : [column row SLC], [2nd dimension, 1st dimension, 3rd dimension]
%   header : the ftk header format, including positionPatient (tx ty tz) and orientationPatient (row vector; col vector; slc vector)
%   figh = figure handle
%   transparency = alpha for the image slice. 
%                   Use values < 1 to get partial transparency when
%                   plotting multiple slices in one figure
%   Note: all access starts from 0
%======================================================================

% as the dicom points to the pixel center, the half pixel offset is needed
xlim3D(1) = header.spacingX*-0.5;
xlim3D(2) = header.spacingX*(header.sizeX-0.5);

ylim3D(1) = header.spacingY*-0.5;
ylim3D(2) = header.spacingY*(header.sizeY-0.5);

showAxis = 1;
if ( ~ishandle(h) )
    figure;
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0.02 0.05 0.95 0.95], 'Renderer', 'zbuffer');
    h = gca;
    
    axis(h, 'ij', 'equal');
    if ( ~showAxis )
        axis(h, 'off');
        set(h, 'Position', [0.2 0.2 0.7 0.7]);
    else
        set(h, 'Position', [0.2 0.2 0.7 0.7]);
    end

    box(h, 'off');
    set(h, 'View', [45 45]);
    set(h, 'XAxisLocation', 'top', 'YAxisLocation', 'left'); 

%     set(h, 'XLim', xlim3D );
%     set(h, 'YLim', ylim3D );
    set(h, 'XColor', [1 0 0], 'YColor', [0 0 1], 'ZColor', [0 0 0]);    
    colormap(gray)
else
    set(gcf, 'Renderer', 'zbuffer');
    colormap(gray);
end

mag = data; 
lowR = 0;
highR = 4096;
mag = normalizeImage2Range(mag, lowR, highR);
magWindowed = normalizeWindowSetting(mag, centre, width);

% for every 2D slice  
x = [-0.5 header.sizeX-0.5];
y = [-0.5 header.sizeY-0.5];

axes(h);
for slc=1:header.sizeZ

    z = slc-1;
    
    % top left
    [wx, wy, wz ] = Image2WorldMrFtk(header, x(1), y(1), z);
    xlim3D(1,1) = wx;
    ylim3D(1,1) = wy;
    zlim3D(1,1) = wz;
    
    % top right
    [wx, wy, wz ] = Image2WorldMrFtk(header, x(2), y(1), z);
    xlim3D(1,2) = wx;
    ylim3D(1,2) = wy;
    zlim3D(1,2) = wz;
    
    % bottom left
    [wx, wy, wz ] = Image2WorldMrFtk(header, x(1), y(2), z);
    xlim3D(2,1) = wx;
    ylim3D(2,1) = wy;
    zlim3D(2,1) = wz;
    
    % bottom right
    [wx, wy, wz ] = Image2WorldMrFtk(header, x(2), y(2), z);
    xlim3D(2,2) = wx;
    ylim3D(2,2) = wy;
    zlim3D(2,2) = wz;
    
    vi = magWindowed(:,:,slc);
    xyImageHandle = surface('XData',xlim3D,'YData',ylim3D,'ZData',zlim3D,'CData', double(vi),'FaceColor','texturemap','EdgeColor','none');
end

return

