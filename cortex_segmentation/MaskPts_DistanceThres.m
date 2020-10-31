

function pts_masked = MaskPts_DistanceThres(mask, header, pts, Thres_MaskingDistance)
% pts in the image coordinate
% the target == 1 in mask and the background 0
% pts (N*4)

num = size(pts, 1);

inds = find(mask>0);
[i, j, k] = ind2sub(size(mask), inds);

points = [j-1 i-1 k-1];
points_world = image2world_DITK2(points, header);
points_world = double(points_world);
% points_world(:,2) = -points_world(:,2);

[minDist, nearestPoints] = GetNearestPoints_VTK(pts(:,2:4), points_world);
index = find(minDist<=Thres_MaskingDistance);

pts_masked = pts(index(:), :);

return;