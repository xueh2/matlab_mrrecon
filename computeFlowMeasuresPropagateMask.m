function [phases, flow, meanVelocity, maxVelocity, cardiacoutput, mask, moco] = computeFlowMeasuresPropagateMask(sourceMag, sourceMask, mag, phs, meanRR, FOV, keyFrame, venc)
% [phases, flow, meanVelocity, maxVelocity, cardiacoutput, mask, moco] = computeFlowMeasuresPropagateMask(sourceMag, sourceMask, mag, phs, meanRR, FOV, keyFrame, venc)
% compute flow measures by propagating the mask

PHS = size(mag, 3);
phases = [0.5:1:PHS] * meanRR * 1000 / PHS;
[mask, mocoMask] = propagateFlowVesselMask(sourceMag, sourceMask, mag(:,:,keyFrame+1));
[mask, moco] = computeFlowVesselMask(mag, keyFrame, FOV, [], mask);
[flow, meanVelocity, maxVelocity, cardiacoutput] = computeFlow(mag, phs, venc, meanRR, 2048, FOV, mask);
