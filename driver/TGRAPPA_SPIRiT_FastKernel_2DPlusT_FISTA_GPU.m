
function [unwrappedIm, fullkspace, sensitivityMap] = TGRAPPA_SPIRiT_FastKernel_2DPlusT_FISTA_GPU(underSampledKspace, reductionFactor, option)
% This function performs the dynamic image reconstruction using SPIRiT
% strategy with 2D + T regularization, FISTA optimization
% option.KernelSize : kernel size [k_fe k_pe]
% option.thresReg : regularization threshold for matrix inversion, e.g. 1e-4
% option.rawFilterFE : kspace filter along FE direction
% option.rawFilterFEStrength : filter strength for kspace data, FE
% option.rawFilterPE : kspace filter along PE direction
% option.rawFilterPEStrength : filter strength for kspace data, PE
% option.acsFilterFE : acs filter along FE direction
% option.acsFilterFEStrength : filter strength for acs data, FE
% option.acsFilterPE : kspace filter along PE direction
% option.acsFilterPEStrength : filter strength for acs data, PE
% option.spatialSmoothingKernel : size of spatial filter for coil sensitivity estimation
% option.zeroFilledSize : zero filled size, if no zero-filling will be done, make it empty
% option.kspaceCenterFE : the center point index along FE
% option.KLTSen : whether to filter sensitivity using KLT
% option.numOfModesKeptKLTSen : number of KLT modes kept
% option.kspaceFullRef : if set, this kspace will be used for kernel estimation
% option.maxIterSPIRiT : maximal number of SPIRIT iterations for linear LSQR
% option.wavWeightSPIRiT : wavelet weights
% option.showIterSPIRiT : if 1, print info for every spirit iteration, only for pocs solver
% option.stopThresSPIRiT : stop threshold for spirit, only for pocs solver
% option.pocsFlagSPIRiT : if 1, use pocs solver, otherwise, use the conjugate solver
% option.cgSPIRiTlambda : regularization parameter for linear lsqr, not used
% option.TVWeightSPIRiT : weights for total variation weights
% option.dataWeightSPIRiT : if <0, the acquired points cannot be changed; if >0, the acquired points can be changed
% this dataWeightSPIRiT is the weight for data fidelity cost term
% option.objTollSPIRiT : conjugate solver threshold for reduction of objective function 
% option.continuationStep : after this number of iterations, the wavelet weight will be changed
% option.wavThresRatio : after continuationStep iterations, the wavelet weight will be changed by 1/wavThresRatio
% option.unwarpMethod : implementation type of the spirit operator 
% option.Itnlim : the number of conjugate solver will be min(option.Itnlim, option.maxIterSPIRiT)
% option.kSizeEigenVectorCoilSensitivity : kernel size for eigen-vector based coil sensitivity estimation
% option.percentageEigenVectorCoilSensitivity : percentage of energy kept for eigen-vector based coil sensitivity estimation
% option.centre : window centre for pocs display 
% option.width  : window width for pocs display

% ===============================================================

% prepare the MrRecon
s = size(underSampledKspace);
Nfe = s(1);
Npe = s(2);
numOfCoil = s(3);
numOfFrames = s(4);

zeroFilledSize = option.zeroFilledSize;
KLTSen = option.KLTSen;
numOfModesKeptKLTSen = option.numOfModesKeptKLTSen;

% SPIRiT
kernelSize = option.KernelSizeSPIRiT;
maxIterSPIRiT = option.maxIterSPIRiT;
wavWeightSPIRiT = option.wavWeightSPIRiT;
showIterSPIRiT = option.showIterSPIRiT;
stopThresSPIRiT = option.stopThresSPIRiT;
pocsFlagSPIRiT = option.pocsFlagSPIRiT;
cgSPIRiTlambda = option.cgSPIRiTlambda;
TVWeightSPIRiT = option.TVWeightSPIRiT;
dataWeightSPIRiT = option.dataWeightSPIRiT;
objTollSPIRiT = option.objTollSPIRiT;
temporalScalingFactorSPIRiT = option.temporalScalingFactorSPIRiT;
continuationStep = option.continuationStep;
wavThresRatio = option.wavThresRatio;

if ( isfield(option, 'unwarpMethod') )
    unwarpMethod = option.unwarpMethod;
else
    unwarpMethod = 'fft';
end

if ( isfield(option, 'Itnlim') )
    Itnlim = option.Itnlim;
else
    Itnlim = 5;
end

centre = option.centre;
width = option.width;

% ----------------------------
% detect the sampling lines
sampledLineLoc = detectSampledLinesDynamic(underSampledKspace);

% ----------------------------
% detect the valide range for FE sampled data
sampledRangeFE = detectSampledRangeFE(underSampledKspace);

% ----------------------------
% generate the filter if required
% generate the filter if required
rawFilterFE = generateKSpaceFilter(option.rawFilterFE, option.rawFilterFEStrength, Nfe, sampledRangeFE, option.kspaceCenterFE);
rawFilterPE = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, Npe, [0 Npe], Npe/2+1);

acsFilterFE = generateKSpaceFilter(option.acsFilterFE, option.acsFilterFEStrength, Nfe, sampledRangeFE, option.kspaceCenterFE);
acsFilterPE = generateKSpaceFilter(option.acsFilterPE, option.acsFilterPEStrength, Npe, [0 Npe], Npe/2+1);

% ===============================================================

% prepare the reference kspace       
ref = zeros(Nfe, Npe, numOfCoil, numOfFrames);

if ( isempty(option.kspaceFullRef) )
    for f=1:numOfFrames
        ref(:,:,:,f) = getSlidingWindowKSpace(underSampledKspace, f, numOfBlocksForRef, reductionFactor);
        ref(:,:,:,f) = performRawDataFilter(ref(:,:,:,f), acsFilterFE, acsFilterPE);
    end
else
    ref = option.kspaceFullRef;
    ref = performRawDataFilter(ref, acsFilterFE, acsFilterPE);
end

if ( KLTSen )
    % perform the KL filtering on ref
    refUsed = reshape(ref, [Nfe*Npe numOfCoil numOfFrames]);
    [a,V,D] = KL_Eigenimage(refUsed);
    % keep 3 eigenmodes
    V2 = V;
    V2(:,1:end-numOfModesKeptKLTSen) = 0;
    s = size(a);
    b = (reshape(a, [s(1)*s(2),s(3)] ));
    refUsed = (reshape(b*V2', [s(1),s(2),s(3)] ) ) ;
    ref = reshape(refUsed, [Nfe Npe numOfCoil numOfFrames]);
    clear refUsed a b V V2 D;
end

refFull = computeTemporalMean(underSampledKspace, reductionFactor);
refFull = performRawDataFilter(refFull, acsFilterFE, acsFilterPE);   

% ===============================================================

tstart = tic;

disp('SPIRiT starting ... ');    

% allocate outputs
if ( isempty(zeroFilledSize) )
    fullkspace = zeros([Nfe Npe numOfCoil numOfFrames]);
    unwrappedIm = zeros([Nfe Npe numOfCoil numOfFrames]);
    sensitivityMap = zeros([Nfe Npe numOfCoil numOfFrames]);
else
    fullkspace = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfCoil numOfFrames]);
    unwrappedIm = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfCoil numOfFrames]);
    sensitivityMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfCoil numOfFrames]);

    % rawFilterPE = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, zeroFilledSize(2), [0 zeroFilledSize(2)], zeroFilledSize(2)/2+1);
end

%% compute the coil sensivity
tic
complexImageSen = ifft2c(refFull);
if ( ~isempty(zeroFilledSize) ) 
    complexImageSen = Zero_Padding_Resize_NoFiltering(complexImageSen, zeroFilledSize(1), zeroFilledSize(2));
end
% [combinedImSen, senMap] = adaptiveCoilCombination2DForSensitivity_PKMethod(complexImageSen, option.spatialSmoothingKernel);

if ( strcmp(option.csmMethod, 'Walsh') )
    params.b1map.spatial_smoothing = 1;
    params.b1map.coilmap_smoothing_blocksize = option.spatialSmoothingKernel;
    senMap = calculateB1map(complexImageSen, params);
    
    if ( ~isempty(zeroFilledSize) )
        S = calculateB1map(ifft2c(refFull), params);
    end    
end

if ( strcmp(option.csmMethod, 'Jun') )
    kSize = option.kSizeEigenVectorCoilSensitivity;
    opts.choice = 1;
    opts.percentage = option.percentageEigenVectorCoilSensitivity;
    tic
    [senMap, eigD, tau]=ED_eigen_2D(refFull, kSize, size(complexImageSen), opts);
    toc
    
    if ( ~isempty(zeroFilledSize) )
        [S, eigD, tau]=ED_eigen_2D(refFull, kSize, size(refFull), opts);
    end
end

if ( strcmp(option.csmMethod, 'Souheil') )
    kSize = [7 7];
    tic
    senMap = CoilSensitivity_Souheil_Walsh(complexImageSen, kSize, 0);
    toc
    
    if ( ~isempty(zeroFilledSize) )
        S = CoilSensitivity_Souheil_Walsh(ifft2c(refFull), kSize, 0);
    end
end

if ( ~isempty(zeroFilledSize) )
    S = repmat(S, [1 1 1 numOfFrames]);
end

disp(['sensitivity computation time is ' num2str(toc)]);

%% perform the recon
% estimate the possible recon window width
deviceInfo = gpuDevice;
maxMem = deviceInfo.TotalMemory/1024/1024; % in Mbytes

memUsage = 2*Nfe*Npe*numOfCoil*numOfCoil*8/1024/1024;
numOfFrameMax = floor(0.3*maxMem/memUsage); % maximal usage of GPU is 1Gbytes
overlapWidth = 3;

if ( numOfFrameMax < 6 )
    numOfFrameMax = 6;
end

reconTasks = [];
ind = 1;
while ind+overlapWidth<numOfFrames    
    if ( isempty(reconTasks) )
        reconTasks = {ind:ind+numOfFrameMax-1};
        ind = ind+numOfFrameMax-overlapWidth;
    else
        reconTasks = [reconTasks {ind:ind+numOfFrameMax-1}];
        ind = ind+numOfFrameMax-overlapWidth;
    end    
end

if ( numel(reconTasks) > 1 )
    prevTask = reconTasks{end-1};
    lastTask = reconTasks{end};
    if ( lastTask(end) > numOfFrames )
        offset = lastTask(end) - numOfFrames;
        lastTask = lastTask - offset;
        overlap = prevTask(end) - lastTask(1);
        if ( overlap > 0.5*numOfFrameMax )
            lastTask(1) = lastTask(1) + (overlap-floor(0.5*numOfFrameMax));
        end    
        reconTasks{end} = lastTask;
    end
end

lastTask = reconTasks{end};
if ( lastTask(end) ~= numOfFrames )
    lastTask = [lastTask lastTask(end)+1:numOfFrames];
    reconTasks{end} = lastTask;
end

numOfTasks = numel(reconTasks);

for tt=1:numOfTasks
    disp(['task ' num2str(tt) ' : ' num2str(reconTasks{tt})]);
end
disp('-------------------------------------');

% gFactor
% for the average all case, all pe lines are used as acs and fitting
option.acqPELines  = 1:Npe;
option.acqAcsLines = 1:Npe;
Coef = GRAPPA_SrcDstChannels_Kernel_2D(refFull, refFull, option, option.thresReg);
if ( ~isempty(zeroFilledSize) )
    CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, zeroFilledSize(1), zeroFilledSize(2), reductionFactor);
else
    CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, Nfe, Npe, reductionFactor);
end
unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff(senMap, CoeffInIm);
gFactor = sqrt(sum(unmixCoeff.*conj(unmixCoeff), 3));

% estimate the SPIRiT kernels
smom = tic;

kspaceRecon = zeros(size(underSampledKspace));
accumTimes = zeros(numOfFrames, 1);

kernelAll = zeros(kernelSize(1), kernelSize(2), numOfCoil, numOfCoil, numOfFrames);
for f=1:numOfFrames
    disp(['  --  Frame ' num2str(f)]); 

    refKSpace = ref(:,:,:,f);    

    tkernel = tic;
    kernelAll(:,:,:,:,f) = spiritCalibration(refKSpace,kernelSize,option.thresRegSPIRiT);
    disp(['kernel calibration time is ' num2str(toc(tkernel))]);
end

for tt=1:numOfTasks
    
    disp(['Task ' num2str(tt)]); 
    
    frames = reconTasks{tt};
    numOfFrames = numel(frames);
    underSampledKspaceCurr = underSampledKspace(:,:,:, frames(:));
    
    GOPs = cell(numOfFrames,1);
    for f=1:numOfFrames
        disp(['  --  Frame ' num2str(f)]); 
        tkernel = tic;
        kernel = kernelAll(:,:,:,:,frames(f));        
        GOPs{f} = SPIRiT(kernel, unwarpMethod,[Nfe Npe]);
        disp(['kernel computation time is ' num2str(toc(tkernel))]);
    end
    GOP2DT = SPIRiT2DPlusT(GOPs);
    clear GOPs

    currS = S(:,:,:,frames(:));
    
    if ( isempty(option.kspaceFullRef) )        
        initialKSpace = underSampledKspaceCurr;
    else
        initialKSpace = option.kspaceFullRef(:,:,:, frames(:));
    end

    % perform CG optimization
    % for the linear part, recon every frame independently
%     kspaceLinear = zeros(size(underSampledKspaceCurr));
%     for f=1:numOfFrames
%         f
%         t = tic;
%         [kspaceLinear(:,:,:,f), RESVEC] = cgSPIRiT_GPU(double(underSampledKspaceCurr(:,:,:,f)),GOP2DT.SpiritGOPs{f}, maxIterSPIRiT, cgSPIRiTlambda, double(underSampledKspaceCurr(:,:,:,f)));
%         disp(['Linear recon time is ' num2str(toc(t))]);
%     end

    if ( wavWeightSPIRiT > 0 )
%         params = initRegularizedCGParams();
%         params.xfmWeight = wavWeightSPIRiT;
%         % params.tikWeight = cgSPIRiTlambda;
%         params.tikWeight = 0;
%         params.TVWeight = TVWeightSPIRiT;
%         params.dataWeight = dataWeightSPIRiT;
%         params.continuationStep = continuationStep;
%         params.wavWeightRatio = wavThresRatio;
%         params.Itnlim = maxIterSPIRiT;
%         params.objToll = objTollSPIRiT;
%         params.temporalScalingFactor = temporalScalingFactorSPIRiT;
%         params.show = showIterSPIRiT;
%         params.sparseTransform = RedundantHarrDWT3D(1);
% 
%         % [kspaceLinear, RESVEC] = cgSPIRiT2DPlusT(double(underSampledKspace),GOP2DT, params.Itnlim, cgSPIRiTlambda, double(underSampledKspace));
% 
%         if ( dataWeightSPIRiT > 0 )
%             kspace = cgRegularizedSPIRiTWithLineSearch_2DPlusT_WithoutNullSpace(double(underSampledKspaceCurr), GOP2DT, min(maxIterSPIRiT, Itnlim), double(kspaceLinear), params, centre, width);
%         else                
%             kspace = cgRegularizedSPIRiTWithLineSearch_2DPlusT_GPU(double(underSampledKspaceCurr), GOP2DT, min(maxIterSPIRiT, Itnlim), double(kspaceLinear), params, centre, width);
%             % kspace = cgRegularizedSPIRiTWithLineSearch_2DPlusT_CoilSensitivity_GPU(double(underSampledKspaceCurr), GOP2DT, min(maxIterSPIRiT, Itnlim), double(kspaceLinear), currS, params, centre, width);
%         end
        
        params = initRegularizedCGParams();
        params.sparseTransform = RedundantHarrDWT3D(1);
        params.sparseTransform.temporalSoftThreshScaleRatio = option.temporalScalingFactorSPIRiT;
        
        tolerance = 1e-8;
        mMask = (underSampledKspaceCurr==0);
        Afun = @(x) GOP2DT*(mMask.*x);
        ATfun = @(x) mMask.*(GOP2DT'*x);
        
        b = -(GOP2DT*underSampledKspaceCurr);
        
        vargin = {'initialization', initialKSpace, 'maxiteration', option.Itnlim, 'lambda', option.lambda, 'tolerance', tolerance, 'ratioFISTA', option.ratioFISTA, 'gFactor', sqrt(gFactor)};
        
        if ( option.UseCoilSensitivity )
            vargin = [vargin 'senMap' currS];
        end
        
        if ( option.mocoFlag )
            Im = SoS_TemporalArray(initialKSpace);
            r = 2048 / max(Im(:));
            header = CreateFtkHeaderInfo(Im, [1 1 6]);
            params.mocoer = MotionCorrection2DPlusT(option.mocoParams.interpolator, round(numOfFrames/2), header);        
            [res, params.mocoer] = params.mocoer.performMoCo(r*Im, option.mocoParams);
            vargin = [vargin 'mocoer' {params.mocoer}];
        end
        
        [kspace, nIter, cost] = SolveFISTA_TC_SureThres_GPU(Afun, ATfun, b, params.sparseTransform, params.sparseTransform', vargin);
        
    else
        kspace = initialKSpace;
    end
        
    kspaceRecon(:,:,:,frames(:)) = kspaceRecon(:,:,:,frames(:)) + kspace;
    accumTimes(frames(:)) = accumTimes(frames(:)) + 1;
end

numOfFrames = size(underSampledKspace, 4);
for f=1:numOfFrames
    kspaceRecon(:,:,:,f) = kspaceRecon(:,:,:,f)/accumTimes(f);
end
kspace = kspaceRecon;

kspaceFiltered = kspace;
for f=1:numOfFrames
    kspaceFiltered(:,:,:,f) = performRawDataFilter(kspace(:,:,:,f), rawFilterFE, rawFilterPE);
end
kspace = kspaceFiltered;

if ( ~isempty(zeroFilledSize) ) 
    kspace = ifft2c(Zero_Padding_Resize_NoFiltering(fft2c(kspace), zeroFilledSize(1), zeroFilledSize(2)));
end

fullkspace = kspace;
unwrappedIm = ifft2c(kspace);
sensitivityMap = senMap; 

disp(['recon time is ' num2str(toc(smom))]);
disp(['Total recon time is ' num2str(toc(tstart))]);
