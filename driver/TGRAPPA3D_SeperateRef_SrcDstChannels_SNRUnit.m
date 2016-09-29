
function [unwrappedIm, fullkspace, sensitivityMap, gFactor] = TGRAPPA3D_SeperateRef_SrcDstChannels_SNRUnit(underSampledKspace, ref, reductionFactor, option)
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
refFull = mean(ref,4);

% compute the coil sensivity
tic
refFullForSenMap = performRawDataFilter(refFull, option.filterRefCOL, option.filterRefLIN);
refFullForSenMap = zpad(refFullForSenMap, Nfe, Npe, numOfCoil);

complexImageSen = ifft2c(refFullForSenMap);
% [combinedImSen, senMap] = adaptiveCoilCombination2DForSensitivity_PKMethod(complexImageSen, option.spatialSmoothingKernel);   

params.b1map.spatial_smoothing = 1;
params.b1map.coilmap_smoothing_blocksize = option.spatialSmoothingKernel;
senMap = calculateB1map(complexImageSen, params);

if ( ~isempty(zeroFilledSize) ) 
    senMap = Zero_Padding_Resize_NoFiltering(senMap, zeroFilledSize(1), zeroFilledSize(2));
end

disp(['sensitivity computation time is ' num2str(toc)]);
% plotComplexImageArrayMontage(senMap);

%% allocate outputs
% as ACS lines are not put back to kspace, hybrid grappa is used and coil combination is performed
if ( isempty(zeroFilledSize) )
    fullkspace = zeros([Nfe Npe numOfFrames]);
    unwrappedIm = zeros([Nfe Npe numOfFrames]);
    sensitivityMap = zeros([Nfe Npe numOfCoil numOfFrames]);
else
    fullkspace = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
    unwrappedIm = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
    sensitivityMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfCoil numOfFrames]);
end

%% compute the kernels
tic
Coef = GRAPPA_SrcDstChannels_Kernel_2D(refFull, refFull, option, option.thresReg);
if ( ~isempty(zeroFilledSize) )
    CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, zeroFilledSize(1), zeroFilledSize(2), reductionFactor);
else
    CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, Nfe, Npe, reductionFactor);
end
unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff(senMap, CoeffInIm);
disp(['kernel computation time is ' num2str(toc)]);

gFactor = sqrt(sum(unmixCoeff.*conj(unmixCoeff), 3)) / reductionFactor;

for f=1:numOfFrames

    disp(['Frame ' num2str(f)]);tic

    reduced_k_data = underSampledKspace(:,:,:,f);

    im = HTGRAPPA_SrcDstChannels_ImagedomainReconWithUnmixCoeff(reduced_k_data, unmixCoeff, zeroFilledSize);
    imfft = fft2c(im);
    fullkspace(:,:,f) = imfft;
    unwrappedIm(:,:,f) = ifft2c(imfft);
    sensitivityMap(:,:,:,f) = senMap;
    disp(['recon time is ' num2str(toc)]);
end

disp(['Total recon time is ' num2str(toc(tstart))]);
