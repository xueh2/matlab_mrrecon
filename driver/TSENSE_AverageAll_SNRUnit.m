
function [unwrappedIm, fullkspace, sensitivityMap, gFactor, unwrappedIm_PixelWise_oneSenMap, unwrappedIm_PixelWise_twoSenMap] = TSENSE_AverageAll_SNRUnit(underSampledKspace, reductionFactor, option)
% This function performs the dynamic image reconstruction using IcePAT strategy (Average all)
% it is sense without FOV limitation recon, average all reference, adaptive coil combionation
% option.KernelSize : kernel size [k_fe k_pe]
% option.KernelPattern : e.g. [-3 0 3 6]
% option.OutPattern : e.g. [0 1 2] if acquired points are fitted; or, [1 2] if acquired points remain unchanged
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
% option.niters: number of maximal iterations for the cg sense recon
% option.lambda: threshold for the cg sense recon
% option.thresForSecondSenMap: threshold for the second sensitivity map
% option.alpha: reg constant for sense snr reg
% option.senseSubMean: whether to subtract mean
% ===============================================================
% prepare the MrRecon

if ( ~isempty(option.zeroFilledSize) )
    underSampledKspace = zpad(underSampledKspace, option.zeroFilledSize(1), option.zeroFilledSize(2), size(underSampledKspace, 3), size(underSampledKspace, 4));
end

s = size(underSampledKspace);
Nfe = s(1);
Npe = s(2);
numOfCoil = s(3);
if( numel(s) > 3 )
    numOfFrames = s(4);
else
    numOfFrames = 1;
end
zeroFilledSize = option.zeroFilledSize;
kspaceCenterFE = option.kspaceCenterFE;

% ----------------------------
% for the average all case, all pe lines are used as acs and fitting
option.acqPELines  = 1:Npe;
option.acqAcsLines = 1:Npe;

% ----------------------------
% detect the sampling lines
sampledLineLoc = detectSampledLinesDynamic(underSampledKspace);
if ( size(sampledLineLoc, 2) ~= size(underSampledKspace, 4) )
    underSampledKspace = underSampledKspace(:,:,:,1:size(sampledLineLoc, 2));
    numOfFrames = size(sampledLineLoc, 2);
end

% ----------------------------
% detect the valide range for FE sampled data
sampledRangeFE = detectSampledRangeFE(underSampledKspace);

% ===============================================================

% prepare the reference kspace      
refFull = computeTemporalMean(underSampledKspace, reductionFactor);
if ( option.dstChaThres > 0 )
    [refFull, E_0, V, dstCha2] = coilReduction(refFull, option.dstCha, option.dstChaThres);
    
    underSampledKspace = applyEigenVector(underSampledKspace, V);
    underSampledKspace = underSampledKspace(:,:,numOfCoil-dstCha2+1:numOfCoil,:);
else
    dstCha2 = numOfCoil;
end
disp(['After coil reduciton, ' num2str(dstCha2) ' channels are used ... ']);

numOfCoil = size(refFull, 3);
dstCHA = numOfCoil;

refFull = performRawDataFilter(refFull, option.filterRefCOL, option.filterRefLIN);

% ===============================================================

tstart = tic;

disp('IcePAT starting ... ');    

% allocate outputs
if ( isempty(zeroFilledSize) )
    fullkspace = zeros([Nfe Npe dstCHA numOfFrames]);
    unwrappedIm = zeros([Nfe Npe numOfFrames]);
    unwrappedIm_PixelWise_oneSenMap = zeros([Nfe Npe numOfFrames]);
    unwrappedIm_PixelWise_twoSenMap = zeros([Nfe Npe numOfFrames]);
    sensitivityMap = zeros([Nfe Npe dstCHA 2 numOfFrames]);
else
    fullkspace = zeros([zeroFilledSize(1) zeroFilledSize(2) dstCHA numOfFrames]);
    unwrappedIm = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
    unwrappedIm_PixelWise_oneSenMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
    unwrappedIm_PixelWise_twoSenMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
    sensitivityMap = zeros([zeroFilledSize(1) zeroFilledSize(2) dstCHA 2 numOfFrames]);
end

%% compute the coil sensivity
tic

% if ( ~isempty(zeroFilledSize) ) 
%     refFullUsed = zpad(refFull, zeroFilledSize(1), zeroFilledSize(2), numOfCoil);
% end
% complexImageSen = ifft2c(refFull);
% [combinedImSen, senMap] = adaptiveCoilCombination2DForSensitivity_PKMethod(complexImageSen, option.spatialSmoothingKernel);   

% params.b1map.spatial_smoothing = 1;
% params.b1map.coilmap_smoothing_blocksize = option.spatialSmoothingKernel;
% senMap = calculateB1map(complexImageSen, params);

if ( ~isempty(zeroFilledSize) ) 
    refFullUsed = zpad(refFull, zeroFilledSize(1), zeroFilledSize(2), numOfCoil);
    complexImageSen = ifft2c(refFullUsed);
else
    complexImageSen = ifft2c(refFull);
end

% if ( strcmp(option.csmMethod, 'Walsh') )
%     params.b1map.spatial_smoothing = 1;
%     params.b1map.coilmap_smoothing_blocksize = option.spatialSmoothingKernel;
%     senMap = calculateB1map(complexImageSen, params);
% end

% if ( strcmp(option.csmMethod, 'Jun') )
    kSize = option.kSizeEigenVectorCoilSensitivity;
    opts.choice = 2;
    opts.percentage = option.percentageEigenVectorCoilSensitivity;
    opts.thd = option.thresForSecondSenMap;
    opts.thres = 1e-2/2;
    
%     refSize = size(refFull);
%     k=24;
%     k_center= refFull((refSize(1)/2 +1 - k):(refSize(1)/2 +k), (refSize(2)/2 +1 - k):(refSize(2)/2 +k),:);

    tic
    % [senMap, eigD, tau]=ED_eigen_2D(refFull, kSize, size(complexImageSen), opts);
%     [senMap, eigD, tau]=ED_eigen_2D_parallel_fov(refFull, kSize, size(complexImageSen), opts);
    [senMap, eigD] = SpiritCSM2D(refFull, [7 7], size(complexImageSen), 0.005);
    
%     ee = eigD(:,:,1);
%     thresEig = mean(ee(:)) * 2;
    
    thresEig = 0.01;
    
    ee = eigD(:,:,2);
    ind = find(ee(:)>thresEig);
    
    senMap2 = senMap(:,:,:,2);
    for cha=1:size(senMap, 3)
        currSen = senMap2(:,:,cha);
        currSen(ind(:)) = 0;
        senMap2(:,:,cha) = currSen;
    end
    senMap(:,:,:,2) = senMap2;
    
%     [combinedIm, senMap, eigD] = Spirit2DForSensitivity(complexImageSen, [5 5], 0.01);
    toc
% end

% if ( strcmp(option.csmMethod, 'Souheil') )
%     kSize = [7 7];
%     tic
%     senMap = CoilSensitivity_Souheil_Walsh(complexImageSen, kSize, 0);
%     toc
% end

% kSize = [7 7];
% senMap_Souheil = CoilSensitivity_Souheil_Walsh(complexImageSen, kSize, 0);
% 
% [combinedIm, senMap, eigD] = Spirit2DForSensitivity(complexImageSen, [5 5], 0.1);

disp(['sensitivity computation time is ' num2str(toc)]);
% plotComplexImageArrayMontage(senMap);

%% kernel estimation
% tic
% Coef = GRAPPA_SrcDstChannels_Kernel_2D(refFullSrc, refFull, option, option.thresReg);
% if ( ~isempty(zeroFilledSize) )
%     CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, zeroFilledSize(1), zeroFilledSize(2), reductionFactor);
% else
%     CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, Nfe, Npe, reductionFactor);
% end
% unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff(senMap, CoeffInIm);
% disp(['kernel computation time is ' num2str(toc)]);
% 
% gFactor = sqrt(sum(unmixCoeff.*conj(unmixCoeff), 3));
% 
% for f=1:numOfFrames
% 
%     disp(['Frame ' num2str(f)]);tic
% 
%     reduced_k_data = underSampledKspace(:,:,:,f); 
%     option.acqPELines  = sampledLineLoc(:,f);
%     im = HTGRAPPA_SrcDstChannels_ImagedomainReconWithUnmixCoeff(reduced_k_data, unmixCoeff, zeroFilledSize);
% 
%     % im = im ./gFactor;
% 
%     imfft = fft2c(im);
%     imfft = performRawDataFilter(imfft, rawFilterFE, rawFilterPE);
%        
%     fullkspace(:,:,f) = imfft;
%     unwrappedIm(:,:,f) = ifft2c(imfft);
%     sensitivityMap(:,:,:,f) = senMap; 
% 
%     % plotComplexImage(im, [2.44 2.44 10], 250, 500, 2, 1/24);
% 
%     disp(['recon time is ' num2str(toc)]);
% end

%% ----------------------------------------

%% calling the sense
% in tpat, we can subtract mean
meanKspace = computeMeanKSpace(underSampledKspace);

% grappa regularized sense
Coef = GRAPPA_SrcDstChannels_Kernel_2D(meanKspace, meanKspace, option, option.thresReg);
if ( ~isempty(zeroFilledSize) )
    CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, zeroFilledSize(1), zeroFilledSize(2), reductionFactor);
else
    CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, Nfe, Npe, reductionFactor);
end
CoeffInIm = CoeffInIm/reductionFactor;

aliasedMeanIm = ifft2c(meanKspace);
grappaRegularizedMeanIm = zeros(size(aliasedMeanIm));
for dst=1:dstCHA
    grappaRegularizedMeanIm(:,:,dst) = sum(aliasedMeanIm.*CoeffInIm(:,:,:,dst), 3);
end

meanKspace = fft2c(grappaRegularizedMeanIm);
meanIm = SensitivityCoilCombination(grappaRegularizedMeanIm, senMap(:,:,:,1));
meanIm = meanIm / sqrt(reductionFactor);

[unmix, gmap] = ismrm_calculate_sense_unmixing(reductionFactor, senMap(:,:,:,1), abs(meanIm), option.alpha);
[unmix2, gFactor] = ismrm_calculate_senseWithoutFOV_unmixing(reductionFactor, senMap, abs(meanIm), option.alpha);

if ( ~option.senseSubMean )
    meanKspace(:) = 0;
    meanIm(:) = 0;
end

for f=1:numOfFrames

    disp(['Frame ' num2str(f)]);tic

    reduced_k_data = underSampledKspace(:,:,:,f);
    reduced_k_data(:, sampledLineLoc(:,f), :) = reduced_k_data(:, sampledLineLoc(:,f), :) - meanKspace(:, sampledLineLoc(:,f), :);
        
    % use the iterative sense
%     [ImCombined, Im, eigD2, senMap2, x0] = senseWithoutFOV(reduced_k_data, refFull, option.niter, option.lambda, senMap, []);
%     ImCombined = ImCombined + meanIm;
    ImCombined = meanIm;
    
    % use the pixelwise sense
    aliasedIm = sqrt(reductionFactor)*ifft2c(reduced_k_data);
    
    unWarppedIm_oneSenMap = sum(aliasedIm.*unmix, 3);
    unWarppedIm_oneSenMap = unWarppedIm_oneSenMap + meanIm;
    
    unWarppedIm_twoSenMap = sum(aliasedIm.*unmix2(:,:,:,1), 3);
    unWarppedIm_twoSenMap = unWarppedIm_twoSenMap + meanIm;
    
    unWarppedIm2 = sum(aliasedIm.*unmix2(:,:,:,2), 3); % figure; imagescn(abs(cat(3, ImCombined, unWarppedIm_twoSenMap, unWarppedIm_oneSenMap, unWarppedIm2, unWarppedIm_twoSenMap+unWarppedIm2)), [0 40], [1 5]);
    unWarppedIm2 = unWarppedIm2 + meanIm;
    
    kspace = fft2c(ImCombined);
    kspace = performRawDataFilter(kspace, rawFilterFE, rawFilterPE);    
    if ( ~isempty(zeroFilledSize) ) 
        kspace = fft2c(Zero_Padding_Resize_NoFiltering(ifft2c(kspace), zeroFilledSize(1), zeroFilledSize(2)));
    end    
    unwrappedIm(:,:,f) = ifft2c(kspace);

    kspace = fft2c(unWarppedIm_oneSenMap);
    if ( ~isempty(zeroFilledSize) ) 
        kspace = fft2c(Zero_Padding_Resize_NoFiltering(ifft2c(kspace), zeroFilledSize(1), zeroFilledSize(2)));
    end    
    unwrappedIm_PixelWise_oneSenMap(:,:,f) = ifft2c(kspace);

    kspace = fft2c(unWarppedIm_twoSenMap);
    if ( ~isempty(zeroFilledSize) ) 
        kspace = fft2c(Zero_Padding_Resize_NoFiltering(ifft2c(kspace), zeroFilledSize(1), zeroFilledSize(2)));
    end
    unwrappedIm_PixelWise_twoSenMap(:,:,f) = ifft2c(kspace);
    
    sensitivityMap(:,:,:,:,f) = senMap;
    
    fullkspace(:,:,:,f) = fft2c(senMap(:,:,:,1) .* repmat(unwrappedIm_PixelWise_twoSenMap(:,:,f), [1 1 numOfCoil]));

    disp(['recon time is ' num2str(toc)]);
end

disp(['Total recon time is ' num2str(toc(tstart))]);
