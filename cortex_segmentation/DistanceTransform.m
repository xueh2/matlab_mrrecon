
function D = DistanceTransform(ObjectVolume)
% the distance transform of ObjectVolume
% target 1, background 0

D = zeros(size(ObjectVolume));

[ysize, xsize, zsize] = size(ObjectVolume);

% i, j, k: row, col, depth

inds_inside = find(ObjectVolume==1);
[i_inside, j_inside, k_inside] = ind2sub(size(ObjectVolume), inds_inside);
p_inside = [j_inside i_inside k_inside];

inds_outside = find(ObjectVolume==0);
[i_outside, j_outside, k_outside] = ind2sub(size(ObjectVolume), inds_outside);
p_outside = [j_outside i_outside k_outside];

[minDist, nearestPoints] = GetNearestPoints_VTK(p_inside, p_outside);

% num = length(inds_inside);
D(inds_inside) = minDist;
return;