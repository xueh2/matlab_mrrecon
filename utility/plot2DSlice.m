function h = plot2DSlice(data, header, currentAxes, previousHandle)
% plot the 2D slice in the axes

% point 1
p1 = [0 0 0];
[p1(1), p1(2), p1(3) ] = Image2WorldMrFtk(header, -0.5, -0.5, 0);

% point 2
p2 = [0 0 0];
[p2(1), p2(2), p2(3) ] = Image2WorldMrFtk(header, header.sizeX-0.5, -0.5, 0);

% point 3
p3 = [0 0 0];
[p3(1), p3(2), p3(3) ] = Image2WorldMrFtk(header, header.sizeX-0.5, header.sizeY-0.5, 0);

% point 4
p4 = [0 0 0];
[p4(1), p4(2), p4(3) ] = Image2WorldMrFtk(header, -0.5, header.sizeY-0.5, 0);

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
