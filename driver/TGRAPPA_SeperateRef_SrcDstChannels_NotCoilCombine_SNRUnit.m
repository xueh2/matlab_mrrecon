
function [unwrappedIm, fullkspace, sensitivityMap, gFactor, E_0, V] = TGRAPPA_SeperateRef_SrcDstChannels_NotCoilCombine_SNRUnit(underSampledKspace, ref, reductionFactor, option)
% This function performs the dynamic image reconstruction with seperate reference
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

s = size(underSampledKspace);
Nfe = s(1);
Npe = s(2);
numOfCoil = s(3);
numOfFrames = s(4);
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

% ----------------------------
% generate the filter if required
rawFilterFE = generateKSpaceFilter(option.rawFilterFE, option.rawFilterFEStrength, Nfe, sampledRangeFE, kspaceCenterFE);
rawFilterPE = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, Npe, [0 Npe], Npe/2+1);

acsFilterFE = generateKSpaceFilter(option.acsFilterFE, option.acsFilterFEStrength, Nfe_ref, sampledRangeFE_ref, Nfe_ref/2+1);
acsFilterPE = generateKSpaceFilter(option.acsFilterPE, option.acsFilterPEStrength, Npe_ref, [0 Npe_ref], Npe_ref/2+1);

tstart = tic;
%% prepare the reference kspace      
% if more than one ref kspace is available, compute the mean
refFull = performRawDataFilter(mean(ref,4), acsFilterFE, acsFilterPE);
refFullSrc = refFull;

if ( option.dstChaThres >= 0 )
    [refFull, E_0, V, dstCha2] = coilReduction(refFull, option.dstCha, option.dstChaThres);
else
    E_0 = [];
    V = [];
    dstCha2 = numOfCoil;
end
disp(['After coil reduciton, ' num2str(dstCha2) ' channels are used ... ']);
numOfCoil = size(refFull, 3);

% compute the coil sensivity
tic
complexImageSen = ifft2c(refFull);
if ( ~isempty(zeroFilledSize) ) 
    complexImageSen = Zero_Padding_Resize_NoFiltering(complexImageSen, zeroFilledSize(1), zeroFilledSize(2));
else
    complexImageSen = Zero_Padding_Resize_NoFiltering(complexImageSen, Nfe, Npe);
end
% [combinedImSen, senMap] = adaptiveCoilCombination2DForSensitivity_PKMethod(complexImageSen, option.spatialSmoothingKernel);   

if ( strcmp(option.csmMethod, 'Walsh') )
    params.b1map.spatial_smoothing = 1;
    params.b1map.coilmap_smoothing_blocksize = option.spatialSmoothingKernel;
    senMap = calculateB1map(complexImageSen, params);
end

if ( strcmp(option.csmMethod, 'Jun') )
    kSize = option.kSizeEigenVectorCoilSensitivity;
    opts.choice = 1;
    opts.percentage = option.percentageEigenVectorCoilSensitivity;
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

%% allocate outputs
% as ACS lines are not put back to kspace, hybrid grappa is used and coil combination is performed
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

%% compute the kernels
tic
Coef = GRAPPA_SrcDstChannels_Kernel_2D(refFullSrc, refFull, option, option.thresReg);
if ( ~isempty(zeroFilledSize) )
    CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, zeroFilledSize(1), zeroFilledSize(2), reductionFactor);
else
    CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, Nfe, Npe, reductionFactor);
end
unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff(senMap, CoeffInIm);
disp(['kernel computation time is ' num2str(toc)]);

gFactor = sqrt(sum(unmixCoeff.*conj(unmixCoeff), 3));

for f=1:numOfFrames

    disp(['Frame ' num2str(f)]);tic

    reduced_k_data = underSampledKspace(:,:,:,f);
    
    option.acqPELines = sampledLineLoc(:,f);
    fkspace = GRAPPA_SrcDstChannels_Recon_2D(reduced_k_data, Coef, option);
    fkspace = performRawDataFilter(fkspace, rawFilterFE, rawFilterPE);                

    im = ifft2c(fkspace);
    if ( ~isempty(zeroFilledSize) )
        im = Zero_Padding_Resize_NoFiltering(im, zeroFilledSize(1), zeroFilledSize(2));
        fkspace = fft2c(im);
    end
    
    fullkspace(:,:,:,f) = fkspace;
    unwrappedIm(:,:,:,f) = im;
    sensitivityMap(:,:,:,f) = senMap;
    disp(['recon time is ' num2str(toc)]);
end

disp(['Total recon time is ' num2str(toc(tstart))]);
