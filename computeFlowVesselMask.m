function [mask, moco] = computeFlowVesselMask(mag, keyFrame, FOV, vesselMaskName, vesselMask)
% mask = computeFlowVesselMask(mag, keyFrame, FOV, vesselMaskName, vesselMask)
% compute flow vessel mask

data  = double(mag);

header = CreateFtkHeaderInfo(data, [1 1 1 1]);
strategy = 'Consecutive';
inverse = 1;
initial = 0;
numOfPre = 0;
iters = [128 128 32];
sigma = 6.0;
neighbor = 2.0;
stepDiv = 3;
moreIterInv = 1;
algo = 'GLCC';
volumePreserving = 0;

[moco, dx, dy, invdx, invdy] = Matlab_PerformTemporalMotionCorrection(data, header, keyFrame, strategy, inverse, initial, numOfPre, iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving);

% warp the vessel mask
if ( ~isempty(vesselMaskName) )
    load(vesselMaskName)

    [J, BW] = roifill(data(:,:,keyFrame), ROI_info_table(1).ROI_x_coordinates, ROI_info_table(1).ROI_y_coordinates);
    mask = repmat(BW, [1 1 size(data, 3)]);
else
    mask = repmat(vesselMask, [1 1 size(data, 3)]);
end
mask = Matlab_PerformWarpingSeries2D(double(mask), header, single(invdx), single(invdy), keyFrame, 'BSpline', 5, 1);
