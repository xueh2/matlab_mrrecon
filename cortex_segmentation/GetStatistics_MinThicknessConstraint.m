
function [mean_S, std_S, points, values] = GetStatistics_MinThicknessConstraint(Internal_SDF, External_SDF, scalarfield, voxelsize, minThickness)
% get the surface points with the thickness larger than the minThickness

% get the surface points using the linear interpolation
level = 0;
[surface_points, num] = GetSurfacePoints(Internal_SDF, level);

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
values = surface_points_values(index(:));
mean_S = mean(values);
std_S = std(values);

return;