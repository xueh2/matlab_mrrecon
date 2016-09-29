
function sampledRangeFE = detectSampledRangeFE(underSampledKspace)
% detect the sampled data range for the FE direction
% sampledRangeFE is [indOfFirstSample indOfLastSample]

partialSum = sum(underSampledKspace, 4);
partialSum = sum(partialSum, 3);
partialSum = sum(partialSum, 2);
ind = find(abs(partialSum)>0);
sampledRangeFE = [min(ind) max(ind)];