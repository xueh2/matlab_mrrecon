function h = plot2DDicomSlice(data, header, currentAxes, previousHandle)
% plot the 2D slice with dicom coordinates in the axes
% h = plot2DDicomSlice(data, header, currentAxes, previousHandle)

h = CreateFtkHeaderInfoFrom(data, []);
h.spacingX = header.PixelSpacing(1);
h.spacingY = header.PixelSpacing(2);
h.spacingZ = header.SliceThickness;
h.positionPatient = header.ImagePositionPatient;
h.orientationPatient(1, :) = header.ImageOrientationPatient(1:3);
h.orientationPatient(2, :) = header.ImageOrientationPatient(4:6);
h.orientationPatient(3, :) = cross(header.ImageOrientationPatient(1:3), header.ImageOrientationPatient(4:6));

% point 1
p1 = [0 0 0];
[p1(1), p1(2), p1(3) ] = Image2WorldMrFtk(h, -0.5, -0.5, 0);

% point 2
p2 = [0 0 0];
[p2(1), p2(2), p2(3) ] = Image2WorldMrFtk(h, h.sizeX-0.5, -0.5, 0);

% point 3
p3 = [0 0 0];
[p3(1), p3(2), p3(3) ] = Image2WorldMrFtk(h, h.sizeX-0.5, h.sizeY-0.5, 0);

% point 4
p4 = [0 0 0];
[p4(1), p4(2), p4(3) ] = Image2WorldMrFtk(h, -0.5, h.sizeY-0.5, 0);

x = [p1(1) p2(1); p4(1) p3(1)];
y = [p1(2) p2(2); p4(2) p3(2)];
z = [p1(3) p2(3); p4(3) p3(3)];

axes(currentAxes);
hold on
       
if ( ishandle(previousHandle) )
    set(previousHandle,'XData',x,'YData',y,'ZData',z,'CData',double(data(:, :, 1)),'FaceColor','texturemap','EdgeColor','none');
    h = previousHandle;
else
    h = surface('XData',x,'YData',y,'ZData',z,'CData',double(data(:, :, 1)),'FaceColor','texturemap','EdgeColor','none');
end

hold off
drawnow
