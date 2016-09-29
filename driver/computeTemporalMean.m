
function kspaceMean = computeTemporalMean(underSampledKspace, reductionFactor)        
% ----------------------------------------------------------
% compute the temporal mean
% underSampledKspace : undersampled kspace [Nfe Npe numOfCoils numOfFrames]
% ----------------------------------------------------------

sampling_location = detectSampledLinesIrregular(underSampledKspace);
% sampling_location = detectSampledLinesDynamic(underSampledKspace);

sampledTimes = zeros([size(underSampledKspace,2) 1]);
numOfFrames = size(underSampledKspace, 4);

for f=1:numOfFrames
    sampledTimes(sampling_location{f}) = sampledTimes(sampling_location{f}) + 1;
end
sampledTimes = repmat(sampledTimes', [size(underSampledKspace, 1) 1]);
sampledTimes(find(sampledTimes==0)) = 1;

numOfCoils = size(underSampledKspace, 3);
kspaceMean = sum(underSampledKspace, 4);
for c=1:numOfCoils
    kspaceMean(:,:,c) = kspaceMean(:,:,c) ./ sampledTimes;
end
