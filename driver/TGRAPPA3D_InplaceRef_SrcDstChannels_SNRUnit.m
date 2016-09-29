
function [unwrappedIm, fullkspace, sensitivityMap, gFactor, E_0, V] = TGRAPPA3D_InplaceRef_SrcDstChannels_SNRUnit(underSampledKspace, ref, reductionFactor, reductionFactor3D, option)
% This function performs the dynamic image reconstruction with inplace reference
% it is hybrid grappa recon, average all reference, adaptive coil combionation
% option.KernelSize : kernel size [k_fe k_pe k_par]
% option.KernelPattern : e.g. [-3 0 3 6]
% option.KernelPatternPAR : e.g. [-3 0 3 6]
% option.OutPattern : e.g. [0 1 2] if acquired points are fitted; or, [1 2] if acquired points remain unchanged
% option.OutPatternPAR : e.g. [0 1 2] if acquired points are fitted; or, [1 2] if acquired points remain unchanged
% option.thresReg : regularization threshold for matrix inversion, e.g. 1e-4
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
% ref : reference scans, following the same line index for the kspace
% ===============================================================

%% prepare the MrRecon

s = size(underSampledKspace);
Nfe = s(1);
Npe = s(2);
numOfCoil = s(3);
Npar = s(4);
zeroFilledSize = option.zeroFilledSize;
kspaceCenterFE = option.kspaceCenterFE;

%% as inplace reference, we need to assemble a full kspace for ref
refSum = sum(ref, 1);
refSum = sum(refSum, 3);
refSum = squeeze(refSum);
refSum = abs(refSum);

sampledRangePE = detectSampledRangeFE(refSum);
sampledRangePAR = detectSampledRangePE(refSum);

oriRef = ref;

useACSLines = 1;

% ref: COL LIN CHA PAR
ref= ref(:, sampledRangePE(1):sampledRangePE(2), :, sampledRangePAR(1):sampledRangePAR(2));

%% start recon
tstart = tic;
disp('Grappa starting ... ');    

%% prepare the reference kspace      
refFull = mean(ref,5);
refFullSrc = refFull;
refDst = ref;

if ( option.dstChaThres > 0 )
% if ( 0 )
    [refFull, E_0, V, dstCha2] = coilReduction(refFull, option.dstCha, option.dstChaThres);
    oriRef = applyEigenVector(oriRef, V);
    oriRef = oriRef(:,:,end-dstCha2+1:end,:);
    
    refDst = applyEigenVector(refDst, V);
    refDst = refDst(:,:,end-dstCha2+1:end,:);
else
    E_0 = [];
    V = [];
    dstCha2 = numOfCoil;
end
disp(['After coil reduciton, ' num2str(dstCha2) ' channels are used ... ']);
numOfCoil = size(refDst, 3);
% ===============================================================

refInKSpaceFull = zeros(Nfe, Npe, numOfCoil, Npar);
refInKSpaceFull(:,sampledRangePE(1):sampledRangePE(2),:,sampledRangePAR(1):sampledRangePAR(2)) = refDst;
refInKSpaceFull = performRawDataFilter3D(refInKSpaceFull, option.filterRefCOL, option.filterRefLIN, option.filterRefPAR);

% compute the coil sensivity
tic

% complexImageSen = ifft2c(refInKSpaceFull);
% if ( ~isempty(zeroFilledSize) ) 
%     complexImageSen = Zero_Padding_Resize_NoFiltering(complexImageSen, zeroFilledSize(1), zeroFilledSize(2));
% end
% [combinedImSen, senMap] = adaptiveCoilCombination2DForSensitivity_PKMethod(complexImageSen, option.spatialSmoothingKernel);   

% params.b1map.spatial_smoothing = 1;
% params.b1map.coilmap_smoothing_blocksize = option.spatialSmoothingKernel;
% senMap = calculateB1map(complexImageSen, params);


complexImageSen = ifft3c_Permute(refInKSpaceFull);

if ( strcmp(option.csmMethod, 'Walsh') )
    params.b1map.spatial_smoothing = 1;
    params.b1map.coilmap_smoothing_blocksize = option.spatialSmoothingKernel;
    sensitivityMap = calculateB1map(complexImageSen, params);
end

if ( strcmp(option.csmMethod, 'Jun') )
    kSize = option.kSizeEigenVectorCoilSensitivity;
    opts.choice = 1;
    opts.percentage = option.percentageEigenVectorCoilSensitivity;
    tic
    [sensitivityMap, eigD, tau]=ED_eigen_2D(refInKSpaceFull, kSize, size(complexImageSen), opts);
    toc
end

if ( strcmp(option.csmMethod, 'Souheil') )
    kSize = [7 7];
    tic
    sensitivityMap = CoilSensitivity_Souheil_Walsh(complexImageSen, kSize, 0);
    toc
end

disp(['sensitivity computation time is ' num2str(toc)]);
% plotComplexImageArrayMontage(sensitivityMap);

% allocate outputs
if ( ~useACSLines )
    fullkspace = zeros([Nfe Npe Npar]);
    unwrappedIm = zeros([Nfe Npe Npar]);
else
    fullkspace = zeros([Nfe Npe numOfCoil Npar]);
    unwrappedIm = zeros([Nfe Npe numOfCoil Npar]);
end

%% compute the kernels
if ( isfield(option, 'rangeUsedFE') & ~isempty(option.rangeUsedFE) )
    refFullSrc = refFullSrc(option.rangeUsedFE(1):option.rangeUsedFE(2),:,:,:);
    refFull = refFull(option.rangeUsedFE(1):option.rangeUsedFE(2),:,:,:);
end

acsSrc = refFullSrc;
acsDst = refFull;

tic
Coef = GRAPPA_SrcDstChannels_Kernel_3D(acsSrc, acsDst, option, option.thresReg);
disp(['Grappa 3D kernel estimation time is ' num2str(toc)]);

tic
kernelIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients_3D(Coef, option, Nfe, Npe, Npar, reductionFactor, reductionFactor3D);
disp(['Convert to image domain kernel time is ' num2str(toc)]);

tic
unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff_3D(sensitivityMap, kernelIm);
disp(['Compute unmixing coefficient time is ' num2str(toc)]);

gFactor = squeeze(sqrt(sum(unmixCoeff.*conj(unmixCoeff), 3)));
gFactor = gFactor / reductionFactor / reductionFactor3D;   

%% recon
if ( ~useACSLines )
    unwrappedIm = HTGRAPPA_SrcDstChannels_ImagedomainReconWithUnmixCoeff_3D(underSampledKspace, unmixCoeff);
    fullkspace = fft3c(unwrappedIm);
else
    aliasedIm = ifft3c_Permute(underSampledKspace);
    aliasedIm = permute(aliasedIm, [1 2 4 3]);
    for cha=1:numOfCoil
        unwrappedIm(:,:,cha,:) = reshape(sum(aliasedIm.*kernelIm(:,:,:,:,cha), 4), [Nfe Npe 1 Npar]);
    end
    
    fullkspace = fft3c_Permute(unwrappedIm);
    fullkspace(:,sampledRangePE(1):sampledRangePE(2),:,sampledRangePAR(1):sampledRangePAR(2)) = oriRef(:,sampledRangePE(1):sampledRangePE(2), :, sampledRangePAR(1):sampledRangePAR(2));
    unwrappedIm = ifft3c_Permute(fullkspace);
end

disp(['Total recon time is ' num2str(toc(tstart))]);
