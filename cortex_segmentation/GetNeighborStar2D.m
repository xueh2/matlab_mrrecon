
function [NeighborStar, NeighborStar_Subs] = GetNeighborStar2D(i, j, shape, offsets)
% get the index of voxels in the neighborstar centered at (i, j)
% neighborNum: 4, 8

num = size(offsets, 1);
NeighborStar_Subs = ones(num, 2);

NeighborStar_Subs(:, 1) = i;
NeighborStar_Subs(:, 2) = j;

NeighborStar_Subs = NeighborStar_Subs + offsets;
NeighborStar = sub2ind(shape, NeighborStar_Subs(:,1), NeighborStar_Subs(:,2));

return;