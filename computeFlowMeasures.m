function [phases, flow, meanVelocity, maxVelocity, cardiacoutput, mask, moco] = computeFlowMeasures(mag, phs, meanRR, FOV, keyFrame, venc, mask_file_name)
% mask = computeFlowVesselMask(mag, keyFrame, FOV, vesselMaskName, vesselMask)
% compute flow measures

PHS = size(mag, 3);
[mask, moco] = computeFlowVesselMask(mag, keyFrame, FOV, mask_file_name);
[flow, meanVelocity, maxVelocity, cardiacoutput] = computeFlow(mag, phs, venc, meanRR, 2048, FOV, mask);

phases = [0.5:1:PHS] * meanRR * 1000 / PHS;
