
function [unwrappedIm, fullkspace, sensitivityMap, kspaceLinear] = TGRAPPA_SPIRiT_FastKernel_2DPlusT_AllinOne_MEX(underSampledKspace, reductionFactor, option)
% This function performs the dynamic image reconstruction using SPIRiT with mex files
% strategy with 2D + T regularization and moco enhanced regularization
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
% sampledLineLoc = detectSampledLinesDynamic(underSampledKspace);

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
    acsFilterPE = generateKSpaceFilter(option.acsFilterPE, option.acsFilterPEStrength, size(ref,2), [0 size(ref,2)], size(ref,2)/2+1);
    ref = performRawDataFilter(ref, acsFilterFE, acsFilterPE);
end

if ( KLTSen )
    % perform the KL filtering on ref
    refUsed = reshape(ref, [Nfe*size(ref,2) numOfCoil numOfFrames]);
    [a,V,D] = KL_Eigenimage(refUsed);
    % keep 3 eigenmodes
    V2 = V;
    V2(:,1:end-numOfModesKeptKLTSen) = 0;
    s = size(a);
    b = (reshape(a, [s(1)*s(2),s(3)] ));
    refUsed = (reshape(b*V2', [s(1),s(2),s(3)] ) ) ;
    ref = reshape(refUsed, [Nfe size(ref,2) numOfCoil numOfFrames]);
    clear refUsed a b V V2 D;
end

refFull = computeTemporalMean(underSampledKspace, reductionFactor);
acsFilterPE = generateKSpaceFilter(option.acsFilterPE, option.acsFilterPEStrength, size(refFull,2), [0 size(refFull,2)], size(refFull,2)/2+1);
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

if ( strcmp(option.csmMethod, 'Walsh') )
    params.b1map.spatial_smoothing = 1;
    params.b1map.coilmap_smoothing_blocksize = option.spatialSmoothingKernel;
    senMap = calculateB1map(complexImageSen, params);
    
    if ( ~isempty(zeroFilledSize) )
        S = calculateB1map(ifft2c(refFull), params);
    else
        S = senMap;
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
    else
        S = senMap;
    end
end

if ( strcmp(option.csmMethod, 'Souheil') )
    kSize = [7 7];
    tic
    [senMap, rho, S] = CoilSensitivity_Souheil_Walsh_PowerMethod(complexImageSen, kSize);
    toc
    
    if ( ~isempty(zeroFilledSize) )
        [S, rho, Stmp] = CoilSensitivity_Souheil_Walsh_PowerMethod(ifft2c(refFull), kSize);
    else
        S = senMap;
    end
end

S = repmat(S, [1 1 1 numOfFrames]);

disp(['sensitivity computation time is ' num2str(toc)]);

%% perform the recon
% estimate the possible recon window width
reconTasks{1} = 1:numOfFrames;
numOfTasks = numel(reconTasks);

for tt=1:numOfTasks
    disp(['task ' num2str(tt) ' : ' num2str(reconTasks{tt})]);
end
disp('-------------------------------------');

% estimate the SPIRiT kernels
smom = tic;

kspaceRecon = zeros(size(underSampledKspace));
accumTimes = zeros(numOfFrames, 1);

G_I_Rep = zeros(Nfe, Npe, numOfCoil, numOfCoil, numOfFrames);

disp(['-----------------------------------']); 
if ( option.multiTaskRecon )
    for rep=1:numOfFrames
        refKSpace = ref(:,:,:,rep);
        tkernel = tic;
        header = CreateFtkHeaderInfo(refKSpace, [1 1 1 1]);
        [kernel2, kernelIm2, conjKernelIm2, G_I, conj_G_I] = Matlab_PerformSPIRiTKernelCalibration2D(single(refKSpace), header, single(refKSpace), header, ... 
            option.KernelSizeSPIRiT, [1 1], [Npe Nfe], option.thresRegSPIRiT);
        disp(['kernel calibration time is ' num2str(toc(tkernel))]);
        
        G_I_Rep(:,:,:,:,rep) = G_I;
    end
else
    refKSpace = mean(ref, 4);
    tkernel = tic;
    header = CreateFtkHeaderInfo(refKSpace, [1 1 1 1]);
    [kernel2, kernelIm2, conjKernelIm2, G_I, conj_G_I] = Matlab_PerformSPIRiTKernelCalibration2D(single(refKSpace),  header, single(refKSpace), header, ... 
        option.KernelSizeSPIRiT, [1 1], [Npe Nfe], option.thresRegSPIRiT);
    disp(['kernel calibration time is ' num2str(toc(tkernel))]);
    G_I_Rep = repmat(G_I, [1 1 1 1 numOfFrames]);
end
clear kernel2 kernelIm2 conjKernelIm2 G_I conj_G_I

disp(['-----------------------------------']); 

for tt=1:numOfTasks
    
    disp(['Task ' num2str(tt)]); 
      
    frames = reconTasks{tt};
    numOfFrames = numel(frames);
    underSampledKspaceCurr = underSampledKspace(:,:,:, frames(:));
    
    currS = single(S(:,:,:,frames(:)));
    
    if ( isempty(option.kspaceFullRef) )        
        initialKSpace = underSampledKspaceCurr;
    else
        initialKSpace = option.kspaceFullRef(:,:,:, frames(:));
    end

    header = CreateFtkHeaderInfo(underSampledKspaceCurr, [1 1 1 1]);
    headerKernel = CreateFtkHeaderInfo(G_I_Rep, [1 1 1 1]);
    
    % perform CG optimization
    if ( option.computeSPIRITLinear )
        kspaceLinear = Matlab_PerformSPIRiTLinearRecon2D(single(underSampledKspaceCurr), header, single(G_I_Rep), headerKernel, single(initialKSpace), ...
            option.performCenteredFFT, option.maxIterSPIRiT, option.stopThresSPIRiT, option.printInfo);
    else
        kspaceLinear = initialKSpace;
    end

    if ( wavWeightSPIRiT > 0 )
        params = initRegularizedCGParams();
        params.xfmWeight = wavWeightSPIRiT;
        % params.tikWeight = cgSPIRiTlambda;
        params.tikWeight = 0;
        params.TVWeight = TVWeightSPIRiT;
        params.dataWeight = dataWeightSPIRiT;
        params.continuationStep = continuationStep;
        params.wavWeightRatio = wavThresRatio;
        params.Itnlim = maxIterSPIRiT;
        params.objToll = objTollSPIRiT;
        params.temporalScalingFactor = temporalScalingFactorSPIRiT;
        params.show = showIterSPIRiT;
        params.sparseTransform = RedundantHarrDWT3D(1);
        
        if ( ~option.performReconWithMoCo )        
            numOfWavLevels = 1;
            % TwoDPlusTRecon = 1;
            performMoCoRecon = 0;

            kspace = Matlab_PerformSPIRiTNonLinearRecon2D(single(underSampledKspaceCurr), header, single(G_I_Rep), ... 
                headerKernel, single(kspaceLinear), single(currS), option.performCenteredFFT, min(maxIterSPIRiT, Itnlim), objTollSPIRiT, wavWeightSPIRiT, ...
                numOfWavLevels, TVWeightSPIRiT, dataWeightSPIRiT, temporalScalingFactorSPIRiT, option.TwoDPlusTRecon, performMoCoRecon, ...
                option.UseCoilSensitivity, option.printInfo);            
        else
            params.mocoParams = option.mocoParams;
                    
            Im = SoS_TemporalArray(initialKSpace);
            r = 2048 / max(Im(:));
            header = CreateFtkHeaderInfo(Im, [1 1 6]);
            params.mocoer = MotionCorrection2DPlusT('BSpline', round(numOfFrames/2), header);        
            [res, params.mocoer] = params.mocoer.performMoCo(r*Im, option.mocoParams);

            initialMoCoIm = params.mocoer*ifft2c(initialKSpace);

%             initialKSpace2 = params.mocoer' * initialMoCoIm;
%             
%             plotComplexImage(ifft2c(initialKSpace(:,:,:,4)));
%             plotComplexImage(initialMoCoIm(:,:,:,4));
%             
%             %% a test zone
%             params.mocoParams.sigma = params.mocoParams.sigma / params.mocoParams.sigmaDivRatio;
%             
%             if ( params.mocoParams.sigma < 16 )
%                 params.mocoParams.sigma = 16
%             end
%             
%             Im = SoS_Image_TemporalArray(initialMoCoIm);
%             r = 2048 / max(Im(:));
%             params.mocoerOri = params.mocoer;
%             [res, params.mocoer] = params.mocoer.performMoCo(r*Im, params.mocoParams);
%             initialMoCoIm2 = params.mocoer * initialMoCoIm;
%             
%             params.mocoerOri = params.mocoerOri.computeConcatenatedDeformationField(params.mocoer.dx, params.mocoer.dy, params.mocoer.dxInv, params.mocoer.dyInv);
%             params.mocoer = params.mocoerOri;            
%             
%             initialKSpace2 = params.mocoer' * initialMoCoIm2;
%             plotComplexImage(initialMoCoIm2(:,:,:,4));
%             plotComplexImage(initialKSpace2(:,:,:,4));
            
            %% perform the recon
            if ( dataWeightSPIRiT > 0 )
                if ( ~option.recoMotionCorrectedImage )
                    if ( ~option.UseCoilSensitivity )
                        kspace = cgRegSPIRiTWithLineSearch_2DPlusT_NoNullSpace_MoCoReg_GPU(single(underSampledKspaceCurr), GOP2DT, min(maxIterSPIRiT, Itnlim), single(kspaceLinear), params, centre, width);
                    else
                        kspace = cgRegSPIRiT_2DPlusT_NoNullSpace_MoCoReg_CoilSensitivity_GPU(single(underSampledKspaceCurr), GOP2DT, min(maxIterSPIRiT, Itnlim), single(kspaceLinear), single(currS), params, centre, width);
                    end
                else
                    kspace = cgRegSPIRiT_2DPlusT_NoNullSpace_MoCoRegIm_CoilSensitivity_GPU(single(underSampledKspaceCurr), GOP2DT, min(maxIterSPIRiT, Itnlim), single(initialMoCoIm), single(currS), params, centre, width);
                end
            else                
                if ( ~option.UseCoilSensitivity )
                    kspace = cgRegularizedSPIRiTWithLineSearch_2DPlusT_MoCoReg_GPU(single(underSampledKspaceCurr), GOP2DT, min(maxIterSPIRiT, Itnlim), single(kspaceLinear), params, centre, width);
                else                
                    kspace = cgRegSPIRiTWithLineSearch_2DPlusT_MoCoReg_CoilSensitivity_GPU(single(underSampledKspaceCurr), GOP2DT, min(maxIterSPIRiT, Itnlim), single(kspaceLinear), single(currS), params, centre, width);
                end
            end        
        end
        
    else
        kspace = kspaceLinear;
    end
        
    kspaceRecon(:,:,:,frames(:)) = kspaceRecon(:,:,:,frames(:)) + kspace;
    accumTimes(frames(:)) = accumTimes(frames(:)) + 1;
end

numOfFrames = size(underSampledKspace, 4);
for f=1:numOfFrames
    kspaceRecon(:,:,:,f) = kspaceRecon(:,:,:,f)/accumTimes(f);
end
kspace = kspaceRecon;

if ( (dataWeightSPIRiT>0) & option.performReconWithMoCo & option.recoMotionCorrectedImage )
    kspace = fft2c(kspace);
end

kspaceFiltered = kspace;
for f=1:numOfFrames
    kspaceFiltered(:,:,:,f) = performRawDataFilter(kspace(:,:,:,f), rawFilterFE, rawFilterPE);
    kspaceLinear(:,:,:,f) = performRawDataFilter(kspaceLinear(:,:,:,f), rawFilterFE, rawFilterPE);
end
kspace = kspaceFiltered;

if ( ~isempty(zeroFilledSize) ) 
    kspace = ifft2c(Zero_Padding_Resize_NoFiltering(fft2c(kspace), zeroFilledSize(1), zeroFilledSize(2)));
    kspaceLinear = ifft2c(Zero_Padding_Resize_NoFiltering(fft2c(kspaceLinear), zeroFilledSize(1), zeroFilledSize(2)));
end

fullkspace = kspace;
unwrappedIm = ifft2c(kspace);
sensitivityMap = senMap; 

disp(['recon time is ' num2str(toc(smom))]);
disp(['Total recon time is ' num2str(toc(tstart))]);
