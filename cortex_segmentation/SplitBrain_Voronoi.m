
function SplittedBrain = SplitBrain_Voronoi(brainmask, centres)
% split brain into different voronoi area
% this function uses the simplest implementation
% more improvement can be achieved by CVT method

index = find(brainmask>0);
[row, col, depth] = ind2sub(size(brainmask), index);
clear index

numofcentres = size(centres, 1);

SplittedBrain = zeros(size(brainmask), 'uint8');

num = size(row, 1);

numsteps = 8;
pointnum_step = round(num/numsteps);

header = 1;
ender = pointnum_step;

for tt = 1:numsteps
    tt
    pointsA = [ col(header:ender) row(header:ender) depth(header:ender)];
    [minDist, nearestPoints] = GetNearestPoints_VTK(pointsA, centres);
%     clear pointsA minDist
    
    for k = header:ender
        currentP = nearestPoints(k-header+1, :);
        label = 0;
        for m = 1:numofcentres
            if ( (currentP(1)==centres(m,1)) & (currentP(2)==centres(m,2)) & (currentP(3)==centres(m,3)) )
                label = m;
                break;
            end
        end
        SplittedBrain(row(k), col(k), depth(k)) = label;
    end
    
    header = ender+1;
    ender = (tt+1)*pointnum_step;
    if ( tt == numsteps-1 )
        ender = num;
    end
end
return