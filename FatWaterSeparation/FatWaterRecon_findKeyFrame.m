function keyFrame = FatWaterRecon_findKeyFrame(Im)
% function keyFrame = FatWaterRecon_findKeyFrame(Im)
%
% this function find the best key frame for each SET
% first, for set = 1 and echo = 1, find the key frame using minimal SSD strategy
% second, for set = 1, registering all REP of set=2 and echo=1 to the key
% frame of first set. The one with minimal deformation is picked as key
% frame
%
% Inputs:
%   Im : the multi-echo images, [Nfe Npe numOfEcho numOfRep numOfSet]
%
% Output:
%   keyFrame : compute key frame for every SET
%
%     ***************************************
%     *  Hui Xue  (hui-xue@siemens.com)     *
%     ***************************************

if ( numel(size(Im)) == 4 ) 
    % only 1 REP
    keyFrame = zeros(size(Im,5), 1);
    return;
end

if ( size(Im, 4) == 1 ) 
    % only 1 REP
    keyFrame = zeros(size(Im,5), 1);
    return;
end

mag = abs(Im);
r = 2048 / max(mag(:));
mag = r * mag;

Nfe = size(Im, 1);
Npe = size(Im, 2);
numOfRep = size(Im, 4);
keyFrame = zeros(size(Im, 5), 1);

% set 1
imSet = squeeze(mag(:,:,1,:,1));

numOfIm = size(imSet, 3);
ssd = zeros(numOfIm, numOfIm);
for ff=1:numOfIm
    for ff2=1:numOfIm
        diffIm = imSet(:,:,ff)-imSet(:,:,ff2);
        ssd(ff, ff2) = norm(diffIm(:));
    end
end

ssd = sort(ssd, 1);
mssd = median(ssd);
[minssd, kf] = min(mssd);
keyFrame(1) = kf - 1;

% set 2
globalMoCoQuality = zeros(numOfRep, 1);

for t=1:numOfRep
    
    data = zeros(Nfe, Npe, 2);
    data(:,:,1) = imSet(:,:,keyFrame(1)+1);
    data(:,:,2) = mag(:,:,1,t,2);
    
    header = CreateFtkHeaderInfo(data, [1 1 1]);    
    [moco, dx, dy, invDx, invDy] = Matlab_PerformTemporalMotionCorrection(double(data), header, 0, 'Direct', 1, 0, 0, [32 32 32], 12.0, 2.0, 3.0, 1, 'GLCC', 0);

    logJacAll = zeros(size(data));
    for f=1:2
        [meanNorm, maxNorm, meanLogJac, maxLogJac, logJac] = analyzeDeformationField2D(dx(:,:,f), dy(:,:,f), header);
        logJacAll(:,:,f) = logJac;
    end

    deformQuality = zeros(header.sizeZ, 1);
    for f=1:header.sizeZ
        x = logJacAll(:,:,f);        
%         option.lts = 1;
%         [res,raw]=fastmcd(x(:), option);        
%         deformQuality(f) = res.center;
        % deformQuality(f) = median(x(:));
        deformQuality(f) = mean(x(:));
    end

    globalMoCoQuality(t) = sum(deformQuality(:));
end

[bestQuality, kf] = min(globalMoCoQuality);
keyFrame(2) = kf-1;

