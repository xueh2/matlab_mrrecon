
function [unwrappedIm, fullkspace, sensitivityMap] = TGRAPPA_SPIRiT_FastKernel_2DPlusT(underSampledKspace, reductionFactor, option)
% This function performs the dynamic image reconstruction using SPIRiT
% strategy with 2D + T regularization
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

    rawFilterPE = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, zeroFilledSize(2), [0 zeroFilledSize(2)], zeroFilledSize(2)/2+1);
end

%% compute the coil sensivity
tic
complexImageSen = ifft2c(refFull);
if ( ~isempty(zeroFilledSize) ) 
    complexImageSen = Zero_Padding_Resize_NoFiltering(complexImageSen, zeroFilledSize(1), zeroFilledSize(2));
end
[combinedImSen, senMap] = adaptiveCoilCombination2DForSensitivity_PKMethod(complexImageSen, option.spatialSmoothingKernel);

% params.b1map.spatial_smoothing = 1;
% params.b1map.coilmap_smoothing_blocksize = option.spatialSmoothingKernel;
% senMap = calculateB1map(complexImageSen, params);
disp(['sensitivity computation time is ' num2str(toc)]);

%% perform the recon
% estimate the SPIRiT kernels

smom = tic;

GOPs = cell(numOfFrames,1);
for f=1:numOfFrames
    disp(['Frame ' num2str(f)]); 
   
    reduced_k_data = underSampledKspace(:,:,:,f);
    refKSpace = ref(:,:,:,f);    
               
    tkernel = tic;
    kernel = spiritCalibration(refKSpace,kernelSize,option.thresRegSPIRiT);
    GOPs{f} = SPIRiT(kernel, unwarpMethod,[Nfe Npe]);
    disp(['kernel computation time is ' num2str(toc(tkernel))]);
end
GOP2DT = SPIRiT2DPlusT(GOPs);
clear GOPs

if ( isempty(option.kspaceFullRef) )        
    initialKSpace = underSampledKspace;
else
    initialKSpace = option.kspaceFullRef;
end
    
% perform CG optimization
% for the linear part, recon every frame independently
kspaceLinear = zeros(size(underSampledKspace));
for f=1:numOfFrames
    t = tic;
    [kspaceLinear(:,:,:,f), RESVEC] = cgSPIRiT(double(underSampledKspace(:,:,:,f)),GOP2DT.SpiritGOPs{f}, maxIterSPIRiT, cgSPIRiTlambda, double(underSampledKspace(:,:,:,f)));
    disp(['Linear recon time is ' num2str(toc(t))]);
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
    params.show = showIterSPIRiT;
    params.sparseTransform = RedundantHarrDWT3D(1);

    % [kspaceLinear, RESVEC] = cgSPIRiT2DPlusT(double(underSampledKspace),GOP2DT, params.Itnlim, cgSPIRiTlambda, double(underSampledKspace));

    if ( dataWeightSPIRiT > 0 )
        kspace = cgRegularizedSPIRiTWithLineSearch_2DPlusT_WithoutNullSpace(double(underSampledKspace), GOP2DT, min(maxIterSPIRiT, Itnlim), double(kspaceLinear), params, centre, width);
    else                
        kspace = cgRegularizedSPIRiTWithLineSearch_2DPlusT(double(underSampledKspace), GOP2DT, min(maxIterSPIRiT, Itnlim), double(kspaceLinear), params, centre, width);
    end
else
    % perform CG optimization
    [kspace, RESVEC] = cgSPIRiT2DPlusT(double(underSampledKspace), GOP2DT, maxIterSPIRiT, cgSPIRiTlambda, double(underSampledKspace));
end
        
if ( ~isempty(zeroFilledSize) ) 
    kspace = ifft2c(Zero_Padding_Resize_NoFiltering(fft2c(kspace), zeroFilledSize(1), zeroFilledSize(2)));
end

kspaceFiltered = kspace;
for f=1:numOfFrames
    kspaceFiltered(:,:,:,f) = performRawDataFilter(kspace(:,:,:,f), rawFilterFE, rawFilterPE);
end

fullkspace = kspaceFiltered;
unwrappedIm = ifft2DForVolume(kspaceFiltered);
sensitivityMap = senMap; 

disp(['recon time is ' num2str(toc(smom))]);
disp(['Total recon time is ' num2str(toc(tstart))]);
