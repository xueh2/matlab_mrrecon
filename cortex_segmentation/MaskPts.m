
function pts_masked = MaskPts(mask, pts)
% pts in the image coordinate
% the target == 1 in mask and the background 0
% pts (N*4)

num = size(pts, 1);
mask = single(abs(mask));
coefficient = col2row( ComputeCoefficients(mask) );

surface_points_values = zeros(num, 1);

for i = 1:num
    x = pts(i, 2);
    y = pts(i, 3);
    z = pts(i, 4);
    surface_points_values(i) = EvaluateSpline2(coefficient, 3, x, y, z);
end

index = find(surface_points_values>0.5);

pts_masked = pts(index(:), :);
return;