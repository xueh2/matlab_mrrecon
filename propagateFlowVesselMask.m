function [mask2, moco] = propagateFlowVesselMask(keyFrame1, mask1, keyFrame2)
% [mask2, moco] = propagateFlowVesselMask(keyFrame1, mask1, keyFrame2)
% moco keyFrame1 to keyFrame2 and propagate the mask1 to mask2

S = size(keyFrame1);
data = zeros(S(1), S(2), 2);
data(:,:,1) = keyFrame1;
data(:,:,2) = keyFrame2;
header = CreateFtkHeaderInfo(data, [1 1 1 1]);
strategy = 'Direct';
inverse = 1;
initial = 0;
numOfPre = 0;
iters = [64 64 64];
sigma = 6.0;
neighbor = 2.0;
stepDiv = 3;
moreIterInv = 1;
algo = 'GLCC';
volumePreserving = 0;

[moco, dx, dy, invdx, invdy] = Matlab_PerformTemporalMotionCorrection(data, header, 0, strategy, inverse, initial, numOfPre, iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving);


% warp the vessel mask
mask = repmat(mask1, [1 1 2]);
mask = Matlab_PerformWarpingSeries2D(double(mask), header, single(invdx), single(invdy), 0, 'BSpline', 5, 1);
mask2 = mask(:,:,2);
