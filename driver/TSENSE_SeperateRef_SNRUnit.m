
function [unwrappedIm, fullkspace, sensitivityMap, gFactor, unwrappedIm_PixelWise_oneSenMap, unwrappedIm_PixelWise_twoSenMap] = TSENSE_SeperateRef_SNRUnit(underSampledKspace, ref, reductionFactor, option)
% This function performs the dynamic image reconstruction using IcePAT strategy (Average all)
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
% ref : reference kspace for kernel estimation
% option.niters: number of maximal iterations for the cg sense recon
% option.lambda: threshold for the cg sense recon
% option.thresForSecondSenMap: threshold for the second sensitivity map
% ===============================================================
% prepare the MrRecon

if ( ~isempty(option.zeroFilledSize) )
    underSampledKspace = zpad(underSampledKspace, option.zeroFilledSize(1), option.zeroFilledSize(2), size(underSampledKspace, 3), size(underSampledKspace, 4));
end

s = size(underSampledKspace);
Nfe = s(1);
Npe = s(2);
numOfCoil = s(3);
if ( numel(s) > 3 )
    numOfFrames = s(4);
else
    numOfFrames = 1;
end

if ( mod(Npe, reductionFactor) ~= 0 )
    Npe = Npe + reductionFactor - mod(Npe, reductionFactor);
    underSampledKspace2 = zeros(Nfe, Npe, numOfCoil, numOfFrames);
    underSampledKspace2(:,1:s(2),:,:) = underSampledKspace;
    underSampledKspace = underSampledKspace2;
end


zeroFilledSize = option.zeroFilledSize;
kspaceCenterFE = option.kspaceCenterFE;

Nfe_ref = size(ref,1);
Npe_ref = size(ref,2);
numOfFrames_ref = size(ref,4);

% ----------------------------
% for the seperate reference case, all pe lines are used as acs and fitting
option.acqPELines  = 1:Npe_ref;
option.acqAcsLines = 1:Npe_ref;

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
sampledRangeFE_ref = detectSampledRangeFE(ref);

tstart = tic;
%% prepare the reference kspace      
% if more than one ref kspace is available, compute the mean
refFull = performRawDataFilter(mean(ref,4), option.filterRefCOL, option.filterRefLIN);

% compute the coil sensivity
tic

% refFullForSenMap = zpad(refFull, Nfe, Npe, numOfCoil);
% complexImageSen = ifft2c(refFullForSenMap);

complexImageSen = Zero_Padding_Resize_NoFiltering(ifft2c(refFull), Nfe, Npe);

kSize = option.kSizeEigenVectorCoilSensitivity;
opts.choice = 2;
opts.percentage = option.percentageEigenVectorCoilSensitivity;
opts.thd = option.thresForSecondSenMap;
opts.thres = option.thres;
% tic
% [senMap, eigD, tau]=ED_eigen_2D_parallel_fov(refFull, kSize, size(complexImageSen), opts);
% toc

tic
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

toc
    
disp(['sensitivity computation time is ' num2str(toc)]);
% plotComplexImageArrayMontage(senMap);

%% allocate outputs
% as ACS lines are not put back to kspace, hybrid grappa is used and coil combination is performed
if ( isempty(zeroFilledSize) )
    fullkspace = zeros([Nfe Npe numOfCoil numOfFrames]);
    unwrappedIm = zeros([Nfe Npe numOfFrames]);
    unwrappedIm_PixelWise_oneSenMap = zeros([Nfe Npe numOfFrames]);
    unwrappedIm_PixelWise_twoSenMap = zeros([Nfe Npe numOfFrames]);
    sensitivityMap = zeros([Nfe Npe numOfCoil 2 numOfFrames]);
else
    fullkspace = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfCoil numOfFrames]);
    unwrappedIm = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
    unwrappedIm_PixelWise_oneSenMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
    unwrappedIm_PixelWise_twoSenMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
    sensitivityMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfCoil 2 numOfFrames]);

    rawFilterPE = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, zeroFilledSize(2), [0 zeroFilledSize(2)], zeroFilledSize(2)/2+1);
end

%% compute the kernels
tic
% Coef = GRAPPA_SrcDstChannels_Kernel_2D(refFull, refFull, option, option.thresReg);
% if ( ~isempty(zeroFilledSize) )
%     CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, zeroFilledSize(1), zeroFilledSize(2), reductionFactor);
% else
%     CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, Nfe, Npe, reductionFactor);
% end
% unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff(senMap, CoeffInIm);
% disp(['kernel computation time is ' num2str(toc)]);
% 
% gFactor = sqrt(sum(unmixCoeff.*conj(unmixCoeff), 3));

snrMap = SensitivityCoilCombination(complexImageSen, senMap(:,:,:,1));
snrMap = abs(snrMap); % figure; imagescn(snrMap);
[unmix, gmap] = ismrm_calculate_sense_unmixing(reductionFactor, senMap(:,:,:,1), snrMap, option.alpha);
[unmix2, gFactor] = ismrm_calculate_senseWithoutFOV_unmixing(reductionFactor, senMap, snrMap, option.alpha);

for f=1:numOfFrames

    disp(['Frame ' num2str(f)]);tic

    reduced_k_data = underSampledKspace(:,:,:,f);

    % [ImCombined, Im, eigD, senMap, x0] = senseWithoutFOV(reduced_k_data, refFull, option.niter, option.lambda, senMap);
    
    % use the pixelwise sense
    aliasedIm = sqrt(reductionFactor)*ifft2c(reduced_k_data);
    
    unWarppedIm_oneSenMap = sum(aliasedIm.*unmix, 3);
    
    unWarppedIm_twoSenMap = sum(aliasedIm.*unmix2(:,:,:,1), 3);  
    unWarppedIm2 = sum(aliasedIm.*unmix2(:,:,:,2), 3);
    
    % figure; imagescn(abs(cat(3, ImCombined, unWarppedIm_twoSenMap, unWarppedIm_oneSenMap, unWarppedIm2, unWarppedIm_twoSenMap+unWarppedIm2)), [0 40], [1 5]);
    
%     kspace = fft2c(ImCombined);
%     kspace = performRawDataFilter(kspace, rawFilterFE, rawFilterPE);    
%     if ( ~isempty(zeroFilledSize) ) 
%         kspace = fft2c(Zero_Padding_Resize_NoFiltering(ifft2c(kspace), zeroFilledSize(1), zeroFilledSize(2)));
%     end    
%     unwrappedIm(:,:,f) = ifft2c(kspace);

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
    unwrappedIm(:,:,f) = ifft2c(kspace);
    
    sensitivityMap(:,:,:,:,f) = senMap;
    
    fullkspace(:,:,:,f) = fft2c(senMap(:,:,:,1) .* repmat(unwrappedIm(:,:,f), [1 1 numOfCoil]));
    
    disp(['recon time is ' num2str(toc)]);
end

disp(['Total recon time is ' num2str(toc(tstart))]);
