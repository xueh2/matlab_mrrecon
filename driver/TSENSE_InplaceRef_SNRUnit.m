
function [unwrappedIm, fullkspace, sensitivityMap, gFactor, unwrappedIm_PixelWise_oneSenMap, unwrappedIm_PixelWise_twoSenMap, E_0, V] = TSENSE_InplaceRef_SNRUnit(underSampledKspace, ref, reductionFactor, option)
% This function performs the dynamic image reconstruction with inplace reference
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
% ref : reference scans, following the same line index for the kspace
% ===============================================================

%% prepare the MrRecon

% if ( ~isempty(option.zeroFilledSize) )
%     underSampledKspace = zpad(underSampledKspace, option.zeroFilledSize(1), option.zeroFilledSize(2), size(underSampledKspace, 3), size(underSampledKspace, 4));
% end

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

%% as inplace reference, we need to assemble a full kspace for ref
sampling_location_ref = detectSampledLines(ref(:,:,:,1));
oriRef = ref;

useACSLines = 1;

if ( sampling_location_ref(2)-sampling_location_ref(1) > 1 )
    ref = ref + underSampledKspace(:,1:size(ref,2),:,:);
end

loc = detectSampledLines(ref(:,:,:,1));
diff = loc([2:end 1]) - loc;
ind = find(diff==1);
startLine = loc(ind(1));
endLine = loc(ind(end))+1;
ref= ref(:,startLine:endLine, :,:);

option.acqAcsLines = 1:size(ref,2);
option.acqPELines = 1:size(ref,2);

%% get the dimensions
Nfe_ref = size(ref,1);
Npe_ref = size(ref,2);

%% detect the sampling lines
sampledLineLoc = detectSampledLinesDynamic(underSampledKspace);
if ( size(sampledLineLoc, 2) ~= size(underSampledKspace, 4) )
    underSampledKspace = underSampledKspace(:,:,:,1:size(sampledLineLoc, 2));
    numOfFrames = size(sampledLineLoc, 2);
end

%% detect the valide range for FE sampled data
sampledRangeFE = detectSampledRangeFE(underSampledKspace);
sampledRangeFE_ref = detectSampledRangeFE(ref);

%% start recon
tstart = tic;

disp('IcePAT starting ... ');    

%% prepare the reference kspace      
% if more than one ref kspace is available, compute the mean
refFull = performRawDataFilter(mean(ref,4), acsFilterFE, acsFilterPE);

% ===============================================================

% prepare the reference kspace      
refFullSrc = refFull;
refDst = ref;

if ( option.dstChaThres > 0 )
% if ( 0 )
    [refFull, E_0, V, dstCha2] = coilReduction(refFull, option.dstCha, option.dstChaThres);
    oriRef = applyEigenVector(oriRef, V);
    oriRef = oriRef(:,:,end-dstCha2+1:end,:);
    
    refDst = applyEigenVector(refDst, V);
    refDst = refDst(:,:,end-dstCha2+1:end,:);
    
    underSampledKspace = applyEigenVector(underSampledKspace, V);
    underSampledKspace = underSampledKspace(:,:,numOfCoil-dstCha2+1:numOfCoil,:);
else
    E_0 = [];
    V = [];
    dstCha2 = numOfCoil;
end
disp(['After coil reduciton, ' num2str(dstCha2) ' channels are used ... ']);
numOfCoil = size(refDst, 3);
% ===============================================================

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

refDst = performRawDataFilter(refDst, option.filterRefCOL, option.filterRefLIN(startLine:endLine));
complexImageSen = ifft2c(refDst);
if ( ~isempty(zeroFilledSize) ) 
    complexImageSen = Zero_Padding_Resize_NoFiltering(complexImageSen, zeroFilledSize(1), zeroFilledSize(2));
else
    complexImageSen = Zero_Padding_Resize_NoFiltering(complexImageSen, Nfe, Npe);
end

kSize = option.kSizeEigenVectorCoilSensitivity;
opts.choice = 2;
opts.percentage = option.percentageEigenVectorCoilSensitivity;
opts.thd = option.thresForSecondSenMap;
opts.thres = option.thres;
tic
[senMap, eigD, tau]=ED_eigen_2D_parallel_fov(refDst, kSize, size(complexImageSen), opts);
toc

disp(['sensitivity computation time is ' num2str(toc)]);
% plotComplexImageArrayMontage(senMap);

% allocate outputs
if ( ~useACSLines )
    % if ACS lines are not put back to kspace, hybrid grappa is used and coil combination is performed
    if ( isempty(zeroFilledSize) )
        fullkspace = zeros([Nfe Npe numOfFrames]);
        unwrappedIm = zeros([Nfe Npe numOfFrames]);
        unwrappedIm_PixelWise_oneSenMap = zeros([Nfe Npe numOfFrames]);
        unwrappedIm_PixelWise_twoSenMap = zeros([Nfe Npe numOfFrames]);
        sensitivityMap = zeros([Nfe Npe numOfCoil numOfFrames]);
    else
        fullkspace = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
        unwrappedIm = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
        unwrappedIm_PixelWise_oneSenMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
        unwrappedIm_PixelWise_twoSenMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
        sensitivityMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfCoil numOfFrames]);

        % rawFilterPE = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, zeroFilledSize(2), [0 zeroFilledSize(2)], zeroFilledSize(2)/2+1);
    end
else
    % if ACS lines are put back to kspace, kspace grappa is used and coil combination is not performed
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

        % rawFilterPE = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, zeroFilledSize(2), [0 zeroFilledSize(2)], zeroFilledSize(2)/2+1);
    end
end

%% compute the kernels

snrMap = SensitivityCoilCombination(complexImageSen, senMap(:,:,:,1));
snrMap = abs(snrMap); % figure; imagescn(snrMap);
[unmix, gmap] = ismrm_calculate_sense_unmixing(reductionFactor, senMap(:,:,:,1), snrMap, option.alpha);
[unmix2, gFactor] = ismrm_calculate_senseWithoutFOV_unmixing(reductionFactor, senMap, snrMap, option.alpha);

for f=1:numOfFrames

    disp(['Frame ' num2str(f)]);tic

    reduced_k_data = underSampledKspace(:,:,:,f);
    reduced_k_data(:,sampling_location_ref(:),:) = oriRef(:,sampling_location_ref(:),:,f);
    
    if ( ~option.useAveAllKernel )
        
        refSrcRep = ref(:,:,:,f);
        refDstRep = refDst(:,:,:,f);
        
        kSize = option.kSizeEigenVectorCoilSensitivity;
        opts.choice = 1;
        opts.percentage = option.percentageEigenVectorCoilSensitivity;
        opts.thd = option.thresForSecondSenMap;
        tic
        [senMap, eigD, tau]=ED_eigen_2D_parallel_fov(refDstRep, kSize, size(complexImageSen), opts);
        toc
    
        disp(['kernel computation time is ' num2str(toc)]);
    end 

    [ImCombined, Im, eigD2, senMap2, x0] = senseWithoutFOV(reduced_k_data, refDst(:,:,:,f), option.niter, option.lambda, senMap);

    % use the pixelwise sense
    aliasedIm = sqrt(reductionFactor)*ifft2c(reduced_k_data);
    
    unWarppedIm_oneSenMap = sum(aliasedIm.*unmix, 3);
    
    unWarppedIm_twoSenMap = sum(aliasedIm.*unmix2(:,:,:,1), 3);       
    unWarppedIm2 = sum(aliasedIm.*unmix2(:,:,:,2), 3);
    
    % figure; imagescn(abs(cat(3, ImCombined, unWarppedIm_twoSenMap, unWarppedIm_oneSenMap, unWarppedIm2, unWarppedIm_twoSenMap+unWarppedIm2)), [0 40], [1 5]);
    
    kspace = fft2c(ImCombined);
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
    
    disp(['recon time is ' num2str(toc)]);
end

disp(['Total recon time is ' num2str(toc(tstart))]);
