
function [unwrappedIm, fullkspace, sensitivityMap] = TGRAPPA_SPIRiT_FastKernel_SrcDst_GPU(underSampledKspace, reductionFactor, option)
% This function performs the dynamic image reconstruction using SrcDst channel spirit and adaptive signal weighting strategy
% -------------------------------------------------------------------------
%                 Parametes for the grappa part
% 
% option.KernelSize : kernel size [k_fe k_pe]
% option.KernelPattern : e.g. [-3 0 3 6]
% option.OutPattern : e.g. [0 1 2] if acquired points are fitted; or, [1 2] if acquired points remain unchanged
% option.thresReg : regularization threshold for matrix inversion, e.g. 1e-4
% option.performGrappa : if true, perform grappa recon
% option.GrappaOnly : if true, only perform grappa recon
%
% -------------------------------------------------------------------------
%                Parameters for the virtual channel recon part
% option.KernelSizeSPIRiT : kernel size [k_fe k_pe] for spirit step
% option.OutKernelSPIRiT : out kernel size of spirit step
% option.thresRegSPIRiT : regularization threshold for matrix inversion, e.g. 1e-4
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
% option.fullKSpaceSrc : full kspace in src channels for kernel estimation
% option.fullKSpaceDst : full kspace in dst channels for kernel estimation

% ------------------------------------------------------------------------
%               Parameters of common
% option.dstCha : the number of kept destination channels
% if option.dstCha < 0, use option.dstChaThres (0.01 e.g.) to determine the number to desination channels
% if option.dstChaThres < 0, all channels are used
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
% option.centre : window centre for pocs display 
% option.width  : window width for pocs display
% =========================================================================

%% prepare the MrRecon

s = size(underSampledKspace);
Nfe = s(1);
Npe = s(2);
numOfCoil = s(3);
numOfFrame = s(4);

zeroFilledSize = option.zeroFilledSize;
KLTSen = option.KLTSen;
numOfModesKeptKLTSen = option.numOfModesKeptKLTSen;

% SPIRiT
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

useGPU = 1;
if ( strcmp(unwarpMethod, 'fft_AllGPU')~=1 )
    useGPU = 0;
end

if ( isfield(option, 'Itnlim') )
    Itnlim = option.Itnlim;
else
    Itnlim = 5;
end

centre = option.centre;
width = option.width;

%% detect the sampling lines
sampledLineLoc = detectSampledLinesDynamic(underSampledKspace);
sampledRangeFE = detectSampledRangeFE(underSampledKspace);

%% perform the grappa recon

if ( option.performGrappa )
    option2 = option;
    option2.dstCha = -1;
    option2.dstChaThres = -1;
    [unwrappedIm_Full_G, fullkspace_Full_G, sensitivityMap_Full_G, gFactor_Full, E_0_Full, V_Full] = TGRAPPA_AverageAll_SrcDstChannels_NotCoilCombine_SNRUnit(underSampledKspace, reductionFactor, option2);
    [unwrappedIm_G, fullkspace_G, sensitivityMap_G, gFactor, E_0, V] = TGRAPPA_AverageAll_SrcDstChannels_NotCoilCombine_SNRUnit(underSampledKspace, reductionFactor, option);

    dstCha = size(fullkspace_G, 3);

    if ( option.GrappaOnly )
        unwrappedIm = unwrappedIm_G;
        fullkspace = fullkspace_G;
        sensitivityMap = sensitivityMap_G;
        return;
    end

    clear unwrappedIm_Full_G unwrappedIm_G sensitivityMap_Full_G sensitivityMap_G gFactor_Full E_0_Full V_Full gFactor option2
else
    fullkspace_Full_G = option.fullKSpaceSrc;
    fullkspace_G = option.fullKSpaceDst;
    E_0 = option.E_0;
    dstCha = size(fullkspace_G, 3);
end

%% generate the filter if required
rawFilterFE = generateKSpaceFilter(option.rawFilterFE, option.rawFilterFEStrength, Nfe, sampledRangeFE, option.kspaceCenterFE);
rawFilterPE = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, Npe, [0 Npe], Npe/2+1);

acsFilterFE = generateKSpaceFilter(option.acsFilterFE, option.acsFilterFEStrength, Nfe, sampledRangeFE, option.kspaceCenterFE);
acsFilterPE = generateKSpaceFilter(option.acsFilterPE, option.acsFilterPEStrength, Npe, [0 Npe], Npe/2+1);

%% for the V-Spirit part, we need to regenerate the undersampled kspace on the virtual channel
underSampledKspaceDst = zeros([Nfe Npe dstCha numOfFrame]);
for f=1:numOfFrame
    underSampledKspaceDst(:,sampledLineLoc(:,f), :, f) = fullkspace_G(:,sampledLineLoc(:,f), :, f);
end

%% perform the normalization
scale_fctr = norm(underSampledKspaceDst(:))/(numOfFrame*sqrt(numOfCoil));
underSampledKspaceDst = underSampledKspaceDst/scale_fctr;
fullkspace_Full_G = fullkspace_Full_G/scale_fctr;
fullkspace_G = fullkspace_G/scale_fctr;
% norm(underSampledKspaceDst(:))
% norm(fullkspace_Full_G(:))
% norm(fullkspace_G(:))

%% if required, perform KL filtering on ref
ref = fullkspace_G;
if ( KLTSen )
    % perform the KL filtering on dst
    refUsed = reshape(fullkspace_G, [Nfe*Npe size(fullkspace_G,3) numOfFrame]);
    [a,V,D] = KL_Eigenimage(refUsed);
    % keep 3 eigenmodes
    V2 = V;
    V2(:,1:end-numOfModesKeptKLTSen) = 0;
    s = size(a);
    b = (reshape(a, [s(1)*s(2),s(3)] ));
    refUsed = (reshape(b*V2', [s(1),s(2),s(3)] ) ) ;
    ref = reshape(refUsed, [Nfe Npe size(fullkspace_G,3) numOfFrame]);
    
    % perform the KL filtering on src
    refUsed = reshape(fullkspace_Full_G, [Nfe*Npe size(fullkspace_Full_G,3) numOfFrame]);
    [a,V,D] = KL_Eigenimage(refUsed);
    % keep 3 eigenmodes
    V2 = V;
    V2(:,1:end-numOfModesKeptKLTSen) = 0;
    s = size(a);
    b = (reshape(a, [s(1)*s(2),s(3)] ));
    refUsed = (reshape(b*V2', [s(1),s(2),s(3)] ) ) ;
    fullkspace_Full_G = reshape(refUsed, [Nfe Npe size(fullkspace_Full_G,3) numOfFrame]);

    clear refUsed a b V V2 D;
end

refFull = mean(fullkspace_G, 4);
refFull = performRawDataFilter(refFull, acsFilterFE, acsFilterPE);

%% ===============================================================

tstart = tic;

disp('V-SPIRiT starting ... ');    

% allocate outputs
if ( isempty(zeroFilledSize) )
    fullkspace = zeros([Nfe Npe dstCha numOfFrame]);
    unwrappedIm = zeros([Nfe Npe dstCha numOfFrame]);
    sensitivityMap = zeros([Nfe Npe dstCha numOfFrame]);
else
    fullkspace = zeros([zeroFilledSize(1) zeroFilledSize(2) dstCha numOfFrame]);
    unwrappedIm = zeros([zeroFilledSize(1) zeroFilledSize(2) dstCha numOfFrame]);
    sensitivityMap = zeros([zeroFilledSize(1) zeroFilledSize(2) dstCha numOfFrame]);

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
for f=1:numOfFrame

    disp(['Frame ' num2str(f)]); smom = tic;
   
    refSrc = fullkspace_Full_G(:,:,:,f);
    reduced_k_data = underSampledKspaceDst(:,:,:,f);
    refDst = ref(:,:,:,f);
    
    optionVSpirit.KernelSize = option.KernelSizeSPIRiT;
    optionVSpirit.OutPattern = option.OutKernelSPIRiT;
    
    tkernel = tic;
    kernelS2D = spiritCalibration_SrcDst_MultipleOutput(refSrc, refDst, optionVSpirit, option.thresRegSPIRiT);
    kernelD2S = spiritCalibration_SrcDst_MultipleOutput(refDst, refSrc, optionVSpirit, option.thresRegSPIRiT);   
    GOP_S2D_D2S = SPIRiTSrcDstCombined(kernelS2D, kernelD2S,unwarpMethod, [Nfe,Npe]);
    disp(['kernel computation time is ' num2str(toc(tkernel))]);

    if pocsFlagSPIRiT 
        % perform POCS optimization
        % not implemented yet
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
            params.reductionFactor = reductionFactor;
            params.gFactor = option.gFactor;            

            if ( ~isempty(E_0) ) 
                params.EigenValue = E_0(end-dstCha+1:end);
                params.virtualChaRatio = sum(E_0(end-dstCha+1:end))/sum(E_0);
            else
                params.EigenValue = [];
                params.virtualChaRatio = 1;
            end
            
            if ( useGPU )
                [kspaceLinear, RESVEC] = cgSPIRiT_SrcDstCombined_GPU(double(reduced_k_data), GOP_S2D_D2S, maxIterSPIRiT, cgSPIRiTlambda, double(reduced_k_data));
                if ( dataWeightSPIRiT > 0 )                    
                    kspace = cgRegSPIRiTWithLineSearch_WithoutNullSpace_SrcDstCombined_GPU(double(reduced_k_data), GOP_S2D_D2S, min(maxIterSPIRiT, Itnlim), double(kspaceLinear), params);
                else                
                    kspace = cgRegSPIRiTWithLineSearch_SrcDstCombined_GPU(double(reduced_k_data), GOP_S2D_D2S, min(maxIterSPIRiT, Itnlim), double(kspaceLinear), params);
                end
            else
                [kspaceLinear, RESVEC] = cgSPIRiT_SrcDstCombined(double(reduced_k_data), GOP_S2D_D2S, maxIterSPIRiT, cgSPIRiTlambda, double(reduced_k_data));
                if ( dataWeightSPIRiT > 0 )                    
                    kspace = cgRegSPIRiTWithLineSearch_WithoutNullSpace_SrcDstCombined(double(reduced_k_data), GOP_S2D_D2S, min(maxIterSPIRiT, Itnlim), double(kspaceLinear), params);
                else                
                    kspace = cgRegSPIRiTWithLineSearch_SrcDstCombined(double(reduced_k_data), GOP_S2D_D2S, min(maxIterSPIRiT, Itnlim), double(kspaceLinear), params);
                end                
            end
        else
            % perform CG optimization
            if ( ~useGPU )
                [kspace, RESVEC] = cgSPIRiT_SrcDstCombined(double(reduced_k_data), GOP_S2D_D2S, maxIterSPIRiT, cgSPIRiTlambda, double(reduced_k_data));
            else
                [kspace, RESVEC] = cgSPIRiT_SrcDstCombined_GPU(double(reduced_k_data), GOP_S2D_D2S, maxIterSPIRiT, cgSPIRiTlambda, double(reduced_k_data));
            end
        end
    end
    
    clear GOP_S2D_D2S
    
    if ( ~isempty(zeroFilledSize) ) 
        kspace = ifft2c(Zero_Padding_Resize_NoFiltering(fft2c(kspace), zeroFilledSize(1), zeroFilledSize(2)));
    end
    kspaceFiltered = performRawDataFilter(kspace, rawFilterFE, rawFilterPE);
    kspaceFiltered = kspaceFiltered*scale_fctr;
    
    fullkspace(:,:,:,f) = kspaceFiltered;
    unwrappedIm(:,:,:,f) = ifft2DForVolume(kspaceFiltered);
    sensitivityMap(:,:,:,f) = senMap; 

    disp(['recon time is ' num2str(toc(smom))]);
end

disp(['Total recon time is ' num2str(toc(tstart))]);
