
function [refKSpace, ref] = estimateRefUsingGeneralMatrixInversion(kspace, kspaceForMoCo, keyFrame, sampledLineLoc, reductionFactor, voxelsize, iters, sigma)
% perform the motion correction to improve the reference estimation

% im = ifft2DForVolume(kspaceForMoCo);
% r = max(abs(im(:)));
% ratio = round(1e4/r)
ratio = 1e5;
imForMoCo = ratio * SoS_TemporalArray(kspaceForMoCo); 
imForMoCo = abs(imForMoCo);
header = CreateFtkHeaderInfo(imForMoCo, voxelsize);
% plotComplexImageArray(imForMoCo, voxelsize, 250, 500, 2, 1/24);

header2D = header;
header2D.sizeZ = 1;
header2D.zsize = 1;

numOfFrame = size(kspace, 4);
numOfFe = size(kspace, 1);
numOfPe = size(kspace, 2);
numOfCoil = size(kspace, 3);

% strategy ='ConsecutiveDeformation';
strategy ='Direct';
inverse = 1;
initial = 0;
numOfPre = 0;
neighbor = 2.0;
stepDiv = 3.0;
moreIterInv = 1;
algo = 'GLCC';
volumePreserving = 0;
suppressSmallDeform = 0;
smallDeformThres = 0.5;

% perform the scout moco to find whether the correction should be performed or not
itersScout = [32 0];
sigmaScout = 128;
[moco, dx, dy, dxInv, dyInv] = Matlab_PerformTemporalMotionCorrection(imForMoCo, header, keyFrame-1, strategy, ...
                                    inverse, initial, numOfPre, itersScout, sigmaScout, neighbor, stepDiv, moreIterInv,...
                                    algo, volumePreserving, suppressSmallDeform, smallDeformThres);
% plotComplexImageArray(moco, voxelsize, 250, 500, 2, 1/24);

meanNorm = zeros(numOfFrame, 1);
maxNorm = zeros(numOfFrame, 1);
meanLogJac = zeros(numOfFrame, 1);
maxLogJac = zeros(numOfFrame, 1);
for f=1:numOfFrame
    [meanNorm(f), maxNorm(f), meanLogJac(f), maxLogJac(f), logJac] = analyzeDeformationField2D(dx(:,:,f), dy(:,:,f), header);
end

disp(['meanNorm ' num2str(max(meanNorm)) ' maxNorm ' num2str(max(maxNorm)) ' meanLogJac ' num2str(max(meanLogJac)) ' maxLogJac ' num2str(max(maxLogJac))]);

if ( (max(maxNorm) < 1.0) & (max(maxLogJac)<0.01) )
    disp('moco ignored ... ');
    refKSpace = sum(kspace, 4); 
    refKSpace = refKSpace / (numOfFrame/reductionFactor);
    ref = ifft2DForVolume(refKSpace);
    return;
end

% -----------------------
% perform the real moco
volumePreserving = 1;
suppressSmallDeform = 1;
smallDeformThres = 0.5;
[moco, dx, dy, dxInv, dyInv] = Matlab_PerformTemporalMotionCorrection(imForMoCo, header, keyFrame-1, strategy, ...
                                    inverse, initial, numOfPre, iters, sigma, neighbor, stepDiv, moreIterInv,...
                                    algo, volumePreserving, suppressSmallDeform, smallDeformThres);

dx = single(dx);
dy = single(dy);
dxInv = single(dxInv);
dyInv = single(dyInv);

% Matlab_SaveAnalyze(single(moco), header, 'moco.hdr');
% Matlab_SaveAnalyze(single(moco-imForMoCo), header, 'diff.hdr');

maxit = 15;
ref = sum(kspace, 4);
ref = ref / (numOfFrame/reductionFactor);
refIm = ifft2DForVolume(ref);
% plotComplexImage(refIm, voxelsize, 250, 500, 2, 1/24)

[Data, E_0, V] = Eigen_Channel(refIm, eye(numOfCoil) );
E_0 = E_0 / sum(E_0);
trivialCoilInd = find(cumsum(E_0)<0.05);
% plot(cumsum(E_0))
lastCoil = trivialCoilInd(end);

sizeIm = size(Data(:,:,1));
for t = 1:numOfFrame
    f = [];
    p = [];
    peInd = sampledLineLoc(:,t);
    nPe = length(peInd);
    for pe=1:nPe
        for fe=1:numOfFe
            f = [f fe];
            p = [p peInd(pe)];
        end
    end

    shots{t} = sub2ind(sizeIm,f,p);
end
    
refMoco = Data;
for c=numOfCoil:-1:lastCoil
    c
    s0 = Data(:,:,c); 
    s = s0;

    tfd.tforms = [];
    tfd.shots  = shots;
    tfd.interpolation = 'bicubic';
    tfd.siz = size(s0);
    for t = 1:numOfFrame
        tfd.deform(t) = struct('dx',dxInv(:,:,t), 'dy',dyInv(:,:,t));
        tfd.deformInv(t) = struct('dx',dx(:,:,t), 'dy',dy(:,:,t));
    end

    tfd.keyShot = keyFrame;
    tfd.header = header2D;
    tfd.nshots = numOfFrame;

    tic
    [sh, flag] = lsqr(@ghosttransmatrixmult_NonRigid,s(:),[],maxit,[],[],[],tfd);
    toc
    sh = reshape(sh,size(s0));  
    refMoco(:,:,c) = sh;    
end
refMoco(:,:,1:lastCoil-1) = Data(:,:,1:lastCoil-1);
refMoco = reshape( refMoco, numOfFe*numOfPe, numOfCoil );
refMoco = refMoco * V';
ref = reshape( refMoco, numOfFe, numOfPe, numOfCoil );
refKSpace = fft2DForVolume(ref);
% plotComplexImage(ref, voxelsize, 250, 500, 2, 1/24)

% ref = zeros(numOfFe, numOfPe, numOfCoil);
% for c=1:numOfCoil
% 
%     p = zeros(size(kspace(:,:,c,1)));
%     for f=1:numOfFrame
%         p(:, sampledLineLoc(:,f)) = p(:, sampledLineLoc(:,f)) + kspace(:, sampledLineLoc(:,f), c, f);
%     end
%     s0 = ifftshift(ifft2(fftshift(p)));
%     % imtool(abs(s0), []);
%     n1 = size(s0,1); n2 = size(s0,2);
%     nshots = numOfFrame;
%     S = zeros(size(s0));
%     for t = 1:nshots
%       f = [];
%       p = [];
%       peInd = sampledLineLoc(:,t);
%       nPe = length(peInd);
%       for pe=1:nPe
%           for fe=1:numOfFe
%               f = [f fe];
%               p = [p peInd(pe)];
%           end
%       end
% 
%       shots{t} = sub2ind(size(s0),f,p);
%     end
% 
%     s = s0;
% 
%     % INPUT FOR ghosttransmatrixmult
%     tfd.tforms = [];
%     tfd.shots  = shots;
%     tfd.interpolation = 'bicubic';
%     tfd.siz = size(s0);
%     for t = 1:nshots
% %         tfd.deform(t) = struct('dx',dx(:,:,t), 'dy',dy(:,:,t));
% %         tfd.deformInv(t) = struct('dx',dxInv(:,:,t), 'dy',dyInv(:,:,t));
%         tfd.deform(t) = struct('dx',dxInv(:,:,t), 'dy',dyInv(:,:,t));
%         tfd.deformInv(t) = struct('dx',dx(:,:,t), 'dy',dy(:,:,t));
%     end
% 
%     tfd.keyShot = keyFrame;
%     tfd.header = header2D;
%     tfd.nshots = nshots;
% 
%     tic
%     [sh, flag] = lsqr(@ghosttransmatrixmult_NonRigid,s(:),[],maxit,[],[],[],tfd);
%     toc
%     sh = reshape(sh,size(s0));
%     % imtool(abs(sh), []);
%     
%     ref(:,:,c) = ifftshift(fft2(fftshift(sh)));
% end

