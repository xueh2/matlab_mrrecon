
function SDF = CreateApproximated_SDF_2D(data)
% create an approximated signed distance function for a 2D bianry volume
% the data should be a volume with isotropic resolution, so that the image coordinates equal to the world coordinate
% the 4-connected neighborhood is used for both inter and outer surfaces
% the distance of surface voxels is assigned to be +/- 0.5 (pixel coordinate)
% the Euler distances of other points are computed using the k-d tree locator
% NOTE: the SDF includes the pixel distance, not in mm

offsets = [-1 0; 1 0; 0 -1; 0 1];

[ysize, xsize] = size(data);

SDF = zeros(size(data), 'double');

points = zeros(xsize*ysize, 3); % store the points on the inter/outer surface

% this loop will get all surface points, so it is safe not to detect all
% boundary points
B = im2col(data,[3 3],'sliding');

index = 0;
for i = 2:xsize-1
    for j = 2:ysize-1
        
        % neighbors = GetNeighbourhood2D(data, i, j, 3);
        neighbors = B(:, (j-1) + (i-2)*(ysize-2) );
        
        hasOne = find(neighbors>0);
        hasZero = find(neighbors==0);

        if ( (length(hasOne)>0) & (length(hasZero)>0) )
           if ( data(j, i) > 0 )
               SDF(j, i) = -0.5;
           else
               SDF(j, i) = 0.5;
           end
           index = index + 1;
           points(index, :) = [i, j, 1];
        end
    end
end
points = points(1:index, :);

% internal points
internalInd = find(data>1);
numOfInternal = length(internalInd);
[j, i] = ind2sub(size(data), internalInd);
points_Internal = [i, j, ones(length(i),1)];
clear internalInd i j

% call the kdtree
% [tmp, tmp, TreeRoot] = kdtree(points, []);
% [ClosestPts, minDist, TreeRoot ] = kdtree([], points_Internal, TreeRoot);
tree = kdtree(points);
[idx, ClosestPts] = kdtree_closestpoint(tree, points_Internal);
minDist = ClosestPts - points_Internal;
minDist = minDist .* minDist;
minDist = sqrt(sum(minDist,2));

% [minDist, nearestPoints] = GetNearestPoints_VTK(points_Internal, points);
for tt = 1:numOfInternal
    i = points_Internal(tt, 1);
    j = points_Internal(tt, 2);
    k = points_Internal(tt, 3);
    SDF(j, i) = -(minDist(tt) + 0.5);
end

% outer points
outerInd = find(data==0);
numOfOuter = length(outerInd);
[j, i] = ind2sub(size(data), outerInd);
points_Outer = [i, j, ones(length(i),1)];
clear outerInd i j

% [tmp, tmp, TreeRoot] = kdtree(points, []);
% [ClosestPts, minDist, TreeRoot ] = kdtree([], points_Outer, TreeRoot);

tree = kdtree(points);
[idx, ClosestPts] = kdtree_closestpoint(tree, points_Outer);
minDist = ClosestPts - points_Outer;
minDist = minDist .* minDist;
minDist = sqrt(sum(minDist,2));

% [minDist, nearestPoints] = GetNearestPoints_VTK(points_Outer, points);
for tt = 1:numOfOuter
    i = points_Outer(tt, 1);
    j = points_Outer(tt, 2);
    k = points_Outer(tt, 3);
    SDF(j, i) = minDist(tt) + 0.5;
end

return