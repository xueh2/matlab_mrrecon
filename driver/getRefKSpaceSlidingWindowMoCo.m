
function [ref, headerRef] = getRefKSpaceSlidingWindowMoCo(underSampledKspace, sampledLineLoc, frame, numOfBlocksForRef, reductionFactor, kspaceFullRef, voxelsize, iters, sigma)
% get the reference kspace for sliding window moco recon
% kspaceFullRef : the full kspace to compute deformation field

Npe = size(underSampledKspace, 2);

currentBlock = floor(frame/reductionFactor) + 1;
if ( mod(frame, reductionFactor) == 0 )
    currentBlock = currentBlock -1;
end

numOfFrameUsed = numOfBlocksForRef*reductionFactor;

if ( mod(numOfFrameUsed, 2) == 0 )
    frameOffsets = -floor(numOfFrameUsed/2):floor(numOfFrameUsed/2)-1;
else
    frameOffsets = -floor(numOfFrameUsed/2):floor(numOfFrameUsed/2);
end
frameInds = frame + frameOffsets;

if( frameInds(1) < 1 )
    frameInds = frameInds - frameInds(1) + 1;
end

if( frameInds(end) > size(underSampledKspace, 4) )
    frameInds = frameInds - frameInds(end) + size(underSampledKspace, 4);
end

ind = frame - frameInds(1) + 1;
disp(['frame ' num2str(frame) ' first ' num2str(frameInds(1)) ' last ' num2str(frameInds(end))]);

% get the kspace
kspace = underSampledKspace(:, :, :, frameInds);
kspaceForMoCo = kspaceFullRef(:,:,:,frameInds);

sampledLineLocBlock = sampledLineLoc(:, frameInds);

[refMoCoKSpace, refMoCo] = estimateRefUsingGeneralMatrixInversion(kspace, kspaceForMoCo, ind, sampledLineLocBlock, reductionFactor, voxelsize, iters, sigma);
ref = refMoCoKSpace;

headerRef = generateDynamicReconHeader(1);
headerRef.sampling_location = 1:Npe;
