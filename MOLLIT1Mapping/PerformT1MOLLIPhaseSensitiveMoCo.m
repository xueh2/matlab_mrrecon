function [im_warpped, dx, dy] = PerformT1MOLLIPhaseSensitiveMoCo(TIs, oriTIs, im, header, isComplexImage)
% if isComplexImage == 1, im is the complex image; otherwise, im is the PSIR image
% im is in the order of oriTIs, not sorted TIs

[TIsorted, ind] = sort(oriTIs);

IRExpNum = zeros(1, numel(oriTIs));
currIRExpNum = 0;
for pp=1:numel(oriTIs)-1
    if ( oriTIs(pp) < oriTIs(pp+1) )
        IRExpNum(pp) = currIRExpNum;
    else
        IRExpNum(pp) = currIRExpNum;
        currIRExpNum = currIRExpNum+1;
    end               
end
    
IRExpNum(end) = currIRExpNum;
IRExpNum = IRExpNum(ind);

AcqTime = 1:numel(oriTIs);
AcqTime = 100 * AcqTime;
AcqTime = AcqTime(ind);
    
% find the corresponding last frame for every IR experiment
lastFrame = zeros(numel(TIs), 1);
for ir=1:numel(TIs)-1
    for l=ir:numel(TIs)

        if ( l==numel(TIs) )
            break;
        end

        if ( oriTIs(l+1) < oriTIs(l) )
            break;
        end
    end
    ind = find(TIsorted==oriTIs(l));
    lastFrame(ir) = ind(1);        
end
lastFrame(end) = numel(TIs);
       
if ( max(abs(im(:))) < 400 )
    r = 2048/max(abs(im(:)));
    im = r*im;
end

[TISorted, ind] = sort(oriTIs);
unwrappedImSorted = zeros(size(im));
for f=1:size(im,3)
    unwrappedImSorted(:,:,f) = im(:,:,ind(f));
end
im = unwrappedImSorted;
           
N = header.sizeZ;

if ( isComplexImage )
    backgroundPhs = im(:,:,end);
    backgroundPhs = backgroundPhs ./ (abs(backgroundPhs)+eps);

    imPS = im .* repmat(conj(backgroundPhs), [1 1 N]);
    imPSReal = real(imPS);
else
    imPSReal = im;
end
      
imMagWithSign = abs(im).*sign(imPSReal);
            
%% perform the moco on the PS image
if ( ~isFileExist('globalMoCoQuality.mat') )
    globalMoCoQuality = zeros(header.sizeZ, 1);
    for t=1:header.sizeZ        
        keyFrame = t-1;
        strategy = 'Consecutive';
        inverse = 1;
        initial = 0;
        numOfPre = 0;
        iters = [32 32 32];
        sigma = 12;
        neighbor = 2.0;
        stepDiv = 3;
        moreIterInv = 1;    
        algo = 'GLCC';
        volumePreserving = 0;
        [moco, dx, dy, invDx, invDy] = Matlab_PerformTemporalMotionCorrection(double(imMagWithSign), header, keyFrame, strategy, inverse, initial, ...
            numOfPre, iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving);

        logJacAll = zeros(size(im));
        for f=1:header.sizeZ
            [meanNorm, maxNorm, meanLogJac, maxLogJac, logJac] = analyzeDeformationField2D(dx(:,:,f), dy(:,:,f), header);
            logJacAll(:,:,f) = logJac;
        end

        deformQuality = zeros(header.sizeZ, 1);
        for f=1:header.sizeZ
            x = logJacAll(:,:,f);        
            option.lts = 1;
            [res,raw]=fastmcd(x(:), option);        
            deformQuality(f) = res.center;
        end

        globalMoCoQuality(t) = sum(deformQuality(:));
    end
    save globalMoCoQuality globalMoCoQuality
else
    load('globalMoCoQuality.mat');
end

[bestQuality, bestKeyFrame] = min(globalMoCoQuality);

keyFrame = bestKeyFrame-1;
strategy = 'Consecutive';
inverse = 1;
initial = 0;
numOfPre = 0;
iters = [32 32 32];
sigma = 24;
neighbor = 2.0;
stepDiv = 3;
moreIterInv = 1;    
algo = 'GLCC';
volumePreserving = 0;
[moco, dx, dy, invDx, invDy] = Matlab_PerformTemporalMotionCorrection(double(imMagWithSign), header, keyFrame, strategy, inverse, initial, ...
    numOfPre, iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving);

meanNormAll = zeros(header.sizeZ,1);
maxNormAll = zeros(header.sizeZ,1);
aveLogJacAll = zeros(header.sizeZ,1);
maxLogJacAll = zeros(header.sizeZ,1);
ind = [];
for f=1:header.sizeZ
    [meanNorm, maxNorm, meanLogJac, maxLogJac, logJac] = analyzeDeformationField2D(dx(:,:,f), dy(:,:,f), header);
    meanNormAll(f) = meanNorm;
    maxNormAll(f) = maxNorm;
    aveLogJacAll(f) = meanLogJac;
    maxLogJacAll(f) = maxLogJac;

    if ( ~((meanLogJac<0.1) & (maxLogJac<0.25) & (meanNorm<5.0) & (maxNorm<8.0)) )
        ind = [ind f];
    end
end

if ( ~isempty(ind) )            
    replacedInd = [];
    for ii=1:numel(ind)               
        if ( ind(ii)==1 | ind(ii)==3 | ind(ii)==2 | ind(ii)==4 )
            dx(:,:,ind(ii)) = dx(:,:,lastFrame(ind(ii)));
            dy(:,:,ind(ii)) = dy(:,:,lastFrame(ind(ii)));
            replacedInd = [replacedInd ind(ii)];
        end
    end

    ind = setdiff(ind, replacedInd); 
end        

% warp the original complex image
if ( isComplexImage )
    im_warpped = PerformComplexImageWarpping(im, header, single(dx), single(dy), keyFrame);
else
    im_warpped = PerformComplexImageWarpping(complex(im, im), header, single(dx), single(dy), keyFrame);
    im_warpped = real(im_warpped);
end

if ( ~isempty(ind) )
    moco(:,:,ind(:)) = imMagWithSign(:,:,ind(:));
    im_warpped(:,:,ind(:)) = im(:,:,ind(:));
end
