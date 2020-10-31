
function [points, num] = RefineSurfacePoints_MinThicknessConstraint_VoxelDistance(pointsInImage, External_SDF, minThickness)
% remove the surface points obtained by vtkPolyData if the thickness is
% below the minThickness
% minThickness: the minimal voxel distance

num = size(pointsInImage, 1);

if ( isempty(External_SDF) == 0 )
    % B-Spline interpolation
    External_SDF = single(abs(External_SDF));
    coefficient = col2row( ComputeCoefficients(External_SDF) );
    surface_points_values = zeros(num, 1);

    for i = 1:num
        x = pointsInImage(i, 1);
        y = pointsInImage(i, 2);
        z = pointsInImage(i, 3);
        surface_points_values(i) = EvaluateSpline2(coefficient, 3, x, y, z);
    end

%     minThickness = minThickness/voxelsize; % mm --> voxel
    index = find(abs(surface_points_values)>=minThickness);

    points = pointsInImage(index(:), :);
else
    points = pointsInImage;
end
num = size(points, 1);

return;