
function [bestOrder_ind, bestOrder, bestNumKnots_ind, bestNumKnots] = PerfSimuData_FindBestOrder_NumberKnots(diffF, OrderBSpline, numOfInternalControlPoints, sigmaStart)

S = size(diffF, 3);

bestOrder = zeros(S, 1);
bestM = zeros(S, 1);

for s=1:S
    diff = diffF(:,:,s);
    [minError, ind] = min(diff(:));
    
    [o, m] = ind2sub(size(diff), ind);
    bestOrder(s) = o;
    bestM(s) = m;   
end

bestNumKnots_ind = median(bestM(sigmaStart(end):end));
bestNumKnots = numOfInternalControlPoints(bestNumKnots_ind);

bestOrder_ind = median(bestOrder(sigmaStart(end):end));
bestOrder = OrderBSpline(bestOrder_ind);
