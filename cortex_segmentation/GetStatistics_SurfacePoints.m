
function [mean_S, std_S, values] = GetStatistics_SurfacePoints(points, scalarfield)
% points are in the image coordinate (in voxel)
% B-Spline interpolation

num = size(points, 1);
scalarfield = single(scalarfield);
coefficient = col2row( ComputeCoefficients(scalarfield) );
values = zeros(num, 1);

for i = 1:num
    x = points(i, 1);
    y = points(i, 2);
    z = points(i, 3);
    values(i) = EvaluateSpline2(coefficient, 3, x, y, z);
end

mean_S = mean(values);
std_S = std(values);

return;