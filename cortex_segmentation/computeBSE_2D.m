
function [BSE, distPercentages, distRealNumbers, dist] = computeBSE_2D(seg2D, manual2D, xvoxelsize, yvoxelsize, distThresholds)
% compute the boundary segmentation errors
% BSE = [mean_BSE std_BSE]
% starting_points, ending_points store the starting and ending coordinates for every short lines
% these points are in pixel coordinates
% distThresholds: number of pixels; count the percentage of distance which is larger than the distThresholds.

sizeSeg = size(seg2D);
sizeManual = size(manual2D);

if ( sum(sizeSeg~=sizeManual)~=0 )
    error('the size of two inputs should be the same ...');
end

num = length(sizeSeg);
size2D = zeros(1,2);
index = 1;
if ( num == 3 )
    for i=1:num
        if ( sizeSeg(i)==1 )
            continue;
        end
        size2D(index) = sizeSeg(i);
        index = index + 1;
    end
    
    seg2D = reshape(seg2D, size2D);
    manual2D = reshape(manual2D, size2D);
end

seg2D(find(seg2D>0)) = 1;
segSDF = CreateApproximated_SDF_2D(seg2D); % the SDF has the pixel coordinates

manual2D(find(manual2D>0)) = 1;
manualSDF = CreateApproximated_SDF_2D(manual2D);

level = 0;
connectivity = 8;

[seg_ps, seg_pe] = CCMS_Contour(segSDF, level, connectivity);
[manual_ps, manual_pe] = CCMS_Contour(manualSDF, level, connectivity);

% compute the symmetric distance using the SDF
num = size(seg_ps, 1);
segDist = zeros(num, 1);

for kk=1:num
    x = seg_ps(kk, 1);
    y = seg_ps(kk, 2);
    value = interp2(manualSDF, x, y);
    segDist(kk) = value;
end

num = size(manual_ps, 1);
manualDist = zeros(num, 1);

for kk=1:num
    x = manual_ps(kk, 1);
    y = manual_ps(kk, 2);
    manualDist(kk) = interp2(segSDF, x, y);
end

dist = [segDist; manualDist];

dist = dist .* (xvoxelsize+yvoxelsize)/2; % pixel, image --> mm, world
numdist = length(dist);

dist = abs(dist);

BSE = [mean(dist) std(dist) max(dist)];

distThresholds = distThresholds .* (xvoxelsize+yvoxelsize)/2;

distPercentages = zeros(size(distThresholds));
distRealNumbers = zeros(size(distThresholds));

absdist = abs(dist);

for i=1:length(distThresholds)
    distRealNumbers(i) = length( find( absdist>distThresholds(i) ) );
    distPercentages(i) = distRealNumbers(i) / numdist;
end

return
