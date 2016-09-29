
function [ref, headerRef] = getRefKSpaceSlidingWindow(underSampledKspace, frame, numOfBlocksForRef, reductionFactor, blockCoord, Nfe, Npe, numOfCoil, minFEUsed, maxFEUsed, sampledLineLoc)        
% ----------------------------------------------------------
% get the reference kspace for sliding window
% underSampledKspace : undersampled kspace [Nfe Npe numOfCoils numOfFrames]
% frame : which frame to compute reference for
% numOfBlocksForRef : number of blocks used to compute the reference
% every reductionFactor consecutive frames is one block
% reductionFactor : regular sampling reductio factor
% ----------------------------------------------------------

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

disp(['frame ' num2str(frame) ' first ' num2str(frameInds(1)) ' last ' num2str(frameInds(end))]);

ref = sum(underSampledKspace(:, :, :, frameInds), 4); 
ref = ref / (numel(frameInds)/reductionFactor);

headerRef = generateDynamicReconHeader(1, reductionFactor, blockCoord, Nfe, Npe, numOfCoil, minFEUsed, maxFEUsed, sampledLineLoc);
headerRef.sampling_location = 1:Npe;
