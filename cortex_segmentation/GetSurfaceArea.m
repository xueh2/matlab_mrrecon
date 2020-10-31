
function [surface_area, cortex_volume] = GetSurfaceArea(Internal, External, level, voxelsize, minThickness, header)
% get surface area and cortex volume
% the minThickness constraint is used
% surface_area: mm^2
% cortex_volume: mL, cm^3

index = find( (Internal>=0) & (External<=0) );
num = length(index);
cortex_volume = num*voxelsize*voxelsize*voxelsize* 10^-3;

[points, num] = GetBoundaryPoints(Internal, 0);

surfaceN = 0;
points2 = zeros(size(points));
for i=1:num
    
    x = points(i,1);
    y = points(i,2);
    z = points(i,3);
    
    TH = External(y, x, z);
    
    if ( abs(TH)>minThickness )
        surfaceN = surfaceN + 1;
        points2(surfaceN,:) = [x y z];
    end
end
points2 = points2(1:surfaceN, :);

surface_area = surfaceN * voxelsize * voxelsize;

return