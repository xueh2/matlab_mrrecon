
function [points, num] = GetSurfacePoints_MinThicknessConstraint(Internal_SDF, External_SDF, level, voxelsize, minThickness)
% get the surface points with the thickness larger than the minThickness
% if the External_SDF is empty matrix, the minThickness constraint is not in use. 

% get the surface points using the linear interpolation
[surface_points, num] = GetSurfacePoints(Internal_SDF, level);

if ( isempty(External_SDF) == 0 )
    % B-Spline interpolation
    External_SDF = single(abs(External_SDF));
    coefficient = col2row( ComputeCoefficients(External_SDF) );
    surface_points_values = zeros(num, 1);

    for i = 1:num
        x = surface_points(i, 1);
        y = surface_points(i, 2);
        z = surface_points(i, 3);
        surface_points_values(i) = EvaluateSpline2(coefficient, 3, x, y, z);
    end

    minThickness = minThickness/voxelsize; % mm --> voxel
    index = find(surface_points_values>=minThickness);

    points = surface_points(index(:), :);
else
    points = surface_points;
end
num = size(points, 1);

return;