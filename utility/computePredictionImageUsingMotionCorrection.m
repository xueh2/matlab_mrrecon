function predictionKSpace = computePredictionImageUsingMotionCorrection(fullkspace, frameInFull, indexOfFrameinFull, iters, sigma)
% compute the prediction complex image array using motion correction

S0 = size(fullkspace);

Nfe = S0(1);
Npe = S0(2);
numOfCoils = S0(3);
numOfFrames = S0(4);

data = zeros(Nfe, Npe, numOfFrames);
data(:,:,1) = SoS(frameInFull);
ind = [1:indexOfFrameinFull-1 indexOfFrameinFull+1:numOfFrames];
data(:,:,2:end) = SoS_TemporalArray(fullkspace(:,:,:,ind));
header = CreateFtkHeaderInfo(data, [1 1 1]);
% plotComplexImageArray(data, voxelsize, newNpe, centre, width, delayTime, 1);

% make the data has sufficient dynamic range
mag = norm(data(:));
data = data * 1e3/mag;

keyFrame = 0;
strategy = 'Direct'
inverse = 1
initial = 0
numOfPre = 0
% iters = [100 100 100]
% sigma = 16
neighbor = 2.0
stepDiv = 3.0
moreIterInv = 1
algo = 'GLCC'
volumePreserving = 0
[moco, dx, dy, invDx, invDy] = PerformTemporalMotionCorrectionComplex(data, header, keyFrame, strategy, inverse, initial, numOfPre, iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving);

% plotComplexImageArray(moco, voxelsize, newNpe, centre, width, delayTime, 1);

% warp the full kspace
invDx = single(invDx);
invDy = single(invDy);
header2D = header;
header2D.sizeZ = 2;

predictionIm = ifft2DForVolume(fullkspace);
predictionIm(:,:,:,indexOfFrameinFull) = ifft2DForVolume(frameInFull);

for f=2:numOfFrames
    
    invDxF = invDx(:,:,f-1:f);
    invDxF(:,:,1) = 0;
    invDyF = invDy(:,:,f-1:f);
    invDyF(:,:,1) = 0;
    
    for c=1:numOfCoils
        realPart = real(predictionIm(:,:,c,[f indexOfFrameinFull]));
        realPart = squeeze(double(realPart));
        moco_real = Matlab_PerformWarpingSeries2D(realPart, header2D, invDxF, invDyF, 0, 'BSpline', 5, 0);

        imagPart = imag(predictionIm(:,:,c,[f indexOfFrameinFull]));
        imagPart = squeeze(double(imagPart));
        moco_imag = Matlab_PerformWarpingSeries2D(imagPart, header2D, invDxF, invDyF, 0, 'BSpline', 5, 0);
        
        if( f <= indexOfFrameinFull )
            predictionIm(:,:,c,f-1) = complex(single(moco_real(:,:,2)), single(moco_imag(:,:,2)));
        else
            predictionIm(:,:,c,f) = complex(single(moco_real(:,:,2)), single(moco_imag(:,:,2)));
        end
    end
end

predictionKSpace = fft2DForVolume(predictionIm);
