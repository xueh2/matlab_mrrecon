
function [unwrappedIm, fullkspace, sensitivityMap] = TGRAPPA_SPIRiT_FastKernel_AverageAll_GPU(underSampledKspace, reductionFactor, option)
% This function performs the dynamic image reconstruction using SPIRiT strategy, averaged-all kernels
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

%% estimate spirit kernel
tkernel = tic;
kernel = spiritCalibration(refFull,kernelSize,option.thresRegSPIRiT);
GOP = SPIRiT(kernel, unwarpMethod,[Nfe Npe]);
disp(['kernel computation time is ' num2str(toc(tkernel))]);

%% perform the recon
for f=1:numOfFrames

    disp(['Frame ' num2str(f)]); smom = tic;
   
    reduced_k_data = underSampledKspace(:,:,:,f);
    refKSpace = ref(:,:,:,f);
    
    if ( isempty(option.kspaceFullRef) )        
        initialKSpace = reduced_k_data;
    else
        initialKSpace = option.kspaceFullRef(:,:,:,f);
    end
               
    if pocsFlagSPIRiT 
        % perform POCS optimization
        kspace = pocsSPIRiT(reduced_k_data, GOP, maxIterSPIRiT, initialKSpace, wavWeightSPIRiT, showIterSPIRiT, stopThresSPIRiT, centre, width, continuationStep, wavThresRatio);
    else
        % perform CG optimization
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
            params.sparseTransform = RedundantHarrDWT2D(2);
            
            [kspaceLinear, RESVEC] = cgSPIRiT_GPU(double(reduced_k_data),GOP, params.Itnlim, cgSPIRiTlambda, double(reduced_k_data));
            
            if ( dataWeightSPIRiT > 0 )
                kspace = cgRegularizedSPIRiTWithLineSearch_WithoutNullSpace_GPU(double(reduced_k_data), GOP, min(maxIterSPIRiT, Itnlim), double(kspaceLinear), params, centre, width);
            else                
                kspace = cgRegularizedSPIRiTWithLineSearch_GPU(double(reduced_k_data), GOP, min(maxIterSPIRiT, Itnlim), double(kspaceLinear), params, centre, width);
            end
        else
            % perform CG optimization
            [kspace, RESVEC] = cgSPIRiT_GPU(double(reduced_k_data),GOP, maxIterSPIRiT, cgSPIRiTlambda, double(reduced_k_data));
        end
    end
        
    if ( ~isempty(zeroFilledSize) ) 
        kspace = ifft2c(Zero_Padding_Resize_NoFiltering(fft2c(kspace), zeroFilledSize(1), zeroFilledSize(2)));
    end
    kspaceFiltered = performRawDataFilter(kspace, rawFilterFE, rawFilterPE);

    fullkspace(:,:,:,f) = kspaceFiltered;
    unwrappedIm(:,:,:,f) = ifft2DForVolume(kspaceFiltered);
    sensitivityMap(:,:,:,f) = senMap; 

    disp(['recon time is ' num2str(toc(smom))]);
end

disp(['Total recon time is ' num2str(toc(tstart))]);
