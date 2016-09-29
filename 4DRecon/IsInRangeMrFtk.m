function inRange = IsInRangeMrFtk(header, wx, wy, wz)
% Find the world point is in the range defined by header or not

[x, y, z ] = World2ImageMrFtk(header, wx, wy, wz);

inRange = 1;

if ( x<0 | x>header.sizeX-1 )
    inRange = 0;
    return;
end

if ( y<0 | y>header.sizeY-1 )
    inRange = 0;
    return;
end

if ( z<0 | z>header.sizeZ-1 )
    inRange = 0;
    return;
end