
function [unwrappedIm, fullkspace, sensitivityMap, gFactor, E_0, V] = TGRAPPA_AverageAll_SrcDstChannels_NotCoilCombine_SNRUnit(underSampledKspace, reductionFactor, option)
% This function performs the dynamic image reconstruction using IcePAT strategy (Average all), but without coil combination
% it is hybrid grappa recon, average all reference, adaptive coil combionation
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
% ===============================================================
% prepare the MrRecon

s = size(underSampledKspace);
Nfe = s(1);
Npe = s(2);
numOfCoil = s(3);
numOfFrames = s(4);
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

% ----------------------------
% generate the filter if required
rawFilterFE = generateKSpaceFilter(option.rawFilterFE, option.rawFilterFEStrength, Nfe, sampledRangeFE, kspaceCenterFE);
rawFilterPE = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, Npe, [0 Npe], Npe/2+1);

acsFilterFE = generateKSpaceFilter(option.acsFilterFE, option.acsFilterFEStrength, Nfe, sampledRangeFE, kspaceCenterFE);
acsFilterPE = generateKSpaceFilter(option.acsFilterPE, option.acsFilterPEStrength, Npe, [0 Npe], Npe/2+1);

% ===============================================================

% prepare the reference kspace      
refFull = computeTemporalMean(underSampledKspace, reductionFactor);
refFull = performRawDataFilter(refFull, acsFilterFE, acsFilterPE);
refFullSrc = refFull;

if ( option.dstCha > 0 )
    [refFull, E_0, V, dstCha2] = coilReduction(refFull, option.dstCha, option.dstChaThres);
elseif ( option.dstChaThres >= 0 )
    [refFull, E_0, V, dstCha2] = coilReduction(refFull, option.dstCha, option.dstChaThres);
else
    E_0 = [];
    V = [];
    dstCha2 = numOfCoil;
end
disp(['After coil reduciton, ' num2str(dstCha2) ' channels are used ... ']);

numOfCoil = size(refFull, 3);
% ===============================================================

tstart = tic;

disp('IcePAT starting ... ');    

% allocate outputs
if ( isempty(zeroFilledSize) )
    fullkspace = zeros([Nfe Npe numOfCoil numOfFrames]);
    unwrappedIm = zeros([Nfe Npe numOfCoil numOfFrames]);
    sensitivityMap = zeros([Nfe Npe numOfCoil numOfFrames]);
else
    fullkspace = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfCoil numOfFrames]);
    unwrappedIm = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfCoil numOfFrames]);
    sensitivityMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfCoil numOfFrames]);
end

% compute the coil sensivity
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
end

if ( strcmp(option.csmMethod, 'Jun') )
    kSize = [5 5];
    opts.choice = 1;
    opts.percentage = 80;
    tic
    [senMap, eigD, tau]=ED_eigen_2D(refFull, kSize, size(complexImageSen), opts);
    toc
end

if ( strcmp(option.csmMethod, 'Souheil') )
    kSize = [7 7];
    tic
    senMap = CoilSensitivity_Souheil_Walsh(complexImageSen, kSize, 0);
    toc
end

disp(['sensitivity computation time is ' num2str(toc)]);
% plotComplexImageArrayMontage(senMap);

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

% for f=1:numOfFrames
% 
%     disp(['Frame ' num2str(f)]);tic
% 
%     reduced_k_data = underSampledKspace(:,:,:,f); 
%     option.acqPELines  = sampledLineLoc(:,f);    
%     kspace = GRAPPA_SrcDstChannels_Recon_2D(reduced_k_data, Coef, option);
% 
%     kspace = performRawDataFilter(kspace, rawFilterFE, rawFilterPE);
%     
%     if ( ~isempty(zeroFilledSize) ) 
%         kspace = fft2c(Zero_Padding_Resize_NoFiltering(ifft2c(kspace), zeroFilledSize(1), zeroFilledSize(2)));
%     end
%     
%     fullkspace(:,:,:,f) = kspace;
%     unwrappedIm(:,:,:,f) = ifft2c(kspace);
%     sensitivityMap(:,:,:,f) = senMap; 
% 
%     disp(['recon time is ' num2str(toc)]);
% end

%% ----------------------------------------

acsSrc = permute(refFullSrc, [2 1 3]);
acsDst = permute(refFull, [2 1 3]);

srcCHA = size(acsSrc, 3);
dstCHA = size(acsDst, 3);

headerSrc = CreateFtkHeaderInfo(acsSrc, [1 1 1 1]);
headerDst = CreateFtkHeaderInfo(acsDst, [1 1 1 1]);

ind = find(option.OutPattern==0);
bFitAcquired = ~isempty(ind);
[kernel, kernelIm, unmixC, gMex] = ... 
    Matlab_PerformVICPAAGKernelCalibrationSrcDst2D(single(acsSrc), headerSrc, single(acsDst), headerDst, ... 
        option.KernelSize(1), option.KernelPattern, 1, 0, [Nfe Npe], option.thresReg, reductionFactor, bFitAcquired);
    
kernelIm = permute(kernelIm, [2 1 3 4]);
gFactor = permute(gMex, [2 1]);

for f=1:numOfFrames

    disp(['Frame ' num2str(f)]);tic

    reduced_k_data = underSampledKspace(:,:,:,f);
    
    aliasedIm = ifft2c(reduced_k_data);
    aliasedIm = repmat(aliasedIm, [1 1 1 dstCHA]);
    unWarppedIm = sum(aliasedIm.*kernelIm, 3);
    unWarppedIm = squeeze(unWarppedIm);
    
    kspace = fft2c(unWarppedIm);
    
    kspace = performRawDataFilter(kspace, rawFilterFE, rawFilterPE);    
    if ( ~isempty(zeroFilledSize) ) 
        kspace = fft2c(Zero_Padding_Resize_NoFiltering(ifft2c(kspace), zeroFilledSize(1), zeroFilledSize(2)));
    end

    fullkspace(:,:,:,f) = kspace;
    unwrappedIm(:,:,:,f) = ifft2c(kspace);
    sensitivityMap(:,:,:,f) = senMap; 

    disp(['recon time is ' num2str(toc)]);
end

disp(['Total recon time is ' num2str(toc(tstart))]);
