
function [ref, headerRef] = getRefKSpaceAverageAll(underSampledKspace, reductionFactor, blockCoord, Nfe, Npe, numOfCoil, minFEUsed, maxFEUsed, sampledLineLoc)        
% ----------------------------------------------------------
% get the reference kspace using average all
% underSampledKspace : undersampled kspace [Nfe Npe numOfCoils numOfFrames]
% ----------------------------------------------------------

Npe = size(underSampledKspace, 2);
ref = computeTemporalMean(underSampledKspace, reductionFactor);
headerRef = generateDynamicReconHeader(1, reductionFactor, blockCoord, Nfe, Npe, numOfCoil, minFEUsed, maxFEUsed, sampledLineLoc);
headerRef.sampling_location = 1:Npe;  
