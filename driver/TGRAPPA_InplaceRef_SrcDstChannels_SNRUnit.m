
function [unwrappedIm, fullkspace, sensitivityMap, gFactor, E_0, V] = TGRAPPA_InplaceRef_SrcDstChannels_SNRUnit(underSampledKspace, ref, reductionFactor, option)
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
refFull = mean(ref,4);

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
else
    E_0 = [];
    V = [];
    dstCha2 = numOfCoil;
end
disp(['After coil reduciton, ' num2str(dstCha2) ' channels are used ... ']);
numOfCoil = size(refDst, 3);
% ===============================================================

refInKSpaceFull = zeros(Nfe, Npe, numOfCoil, numOfFrames);
refInKSpaceFull(:,startLine:endLine,:,:) = refDst;
refInKSpaceFull = performRawDataFilter(refInKSpaceFull, option.filterRefCOL, option.filterRefLIN);

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


complexImageSen = ifft2c(refInKSpaceFull);
if ( ~isempty(zeroFilledSize) ) 
    complexImageSen = Zero_Padding_Resize_NoFiltering(complexImageSen, zeroFilledSize(1), zeroFilledSize(2));
end

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
    [senMap, eigD, tau]=ED_eigen_2D(refInKSpaceFull, kSize, size(complexImageSen), opts);
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

% allocate outputs
if ( ~useACSLines )
    % if ACS lines are not put back to kspace, hybrid grappa is used and coil combination is performed
    if ( isempty(zeroFilledSize) )
        fullkspace = zeros([Nfe Npe numOfFrames]);
        unwrappedIm = zeros([Nfe Npe numOfFrames]);
        sensitivityMap = zeros([Nfe Npe numOfCoil numOfFrames]);
    else
        fullkspace = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
        unwrappedIm = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfFrames]);
        sensitivityMap = zeros([zeroFilledSize(1) zeroFilledSize(2) numOfCoil numOfFrames]);

        % rawFilterPE = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, zeroFilledSize(2), [0 zeroFilledSize(2)], zeroFilledSize(2)/2+1);
    end
else
    % if ACS lines are put back to kspace, kspace grappa is used and coil combination is not performed
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
end

%% compute the kernels

if ( option.useAveAllKernel )
    tic
%     Coef = GRAPPA_SrcDstChannels_Kernel_2D(refFullSrc, refFull, option, option.thresReg);
%     if ( ~isempty(zeroFilledSize) )
%         CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, zeroFilledSize(1), zeroFilledSize(2), reductionFactor);
%     else
%         CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, Nfe, Npe, reductionFactor);
%     end
%     unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff(senMap, CoeffInIm);
%     disp(['kernel computation time is ' num2str(toc)]);
% 
%     gFactor = sqrt(sum(unmixCoeff.*conj(unmixCoeff), 3));
    
    if ( isfield(option, 'rangeUsedFE') & ~isempty(option.rangeUsedFE) )
        refFullSrc = refFullSrc(option.rangeUsedFE(1):option.rangeUsedFE(2),:,:);
        refFull = refFull(option.rangeUsedFE(1):option.rangeUsedFE(2),:,:);
    end

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

    unmixCoeff = permute(unmixC, [2 1 3 4]);
    kernelIm = permute(kernelIm, [2 1 3 4]);
    gFactor = permute(gMex, [2 1]);

end 

R = numel(underSampledKspace(:,:,:,1)) / (numel(find(abs(underSampledKspace(:,:,:,1))>0))+Nfe*numel(sampling_location_ref));
gFactor = gFactor / R;   

for f=1:numOfFrames

    disp(['Frame ' num2str(f)]);tic

    reduced_k_data = underSampledKspace(:,:,:,f);

    if ( ~option.useAveAllKernel )
        
        refSrcRep = ref(:,:,:,f);
        refDstRep = refDst(:,:,:,f);
        
        tic
%         Coef = GRAPPA_SrcDstChannels_Kernel_2D(refSrcRep, refDstRep, option, option.thresReg);
%         if ( ~isempty(zeroFilledSize) )
%             CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, zeroFilledSize(1), zeroFilledSize(2), reductionFactor);
%         else
%             CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, Nfe, Npe, reductionFactor);
%         end
%         unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff(senMap, CoeffInIm);
%         
%         if ( f == 1 ) 
%             gFactor = sqrt(sum(unmixCoeff.*conj(unmixCoeff), 3)); 
%         end
        
        acsSrc = permute(refSrcRep, [2 1 3]);
        acsDst = permute(refDstRep, [2 1 3]);

        srcCHA = size(acsSrc, 3);
        dstCHA = size(acsDst, 3);

        headerSrc = CreateFtkHeaderInfo(acsSrc, [1 1 1 1]);
        headerDst = CreateFtkHeaderInfo(acsDst, [1 1 1 1]);

        ind = find(option.OutPattern==0);
        bFitAcquired = ~isempty(ind);
        [kernel, kernelIm, unmixC, gMex] = ... 
            Matlab_PerformVICPAAGKernelCalibrationSrcDst2D(single(acsSrc), headerSrc, single(acsDst), headerDst, ... 
                option.KernelSize(1), option.KernelPattern, 1, 0, [Nfe Npe], option.thresReg, reductionFactor, bFitAcquired);

        Coef = kernel;
        unmixCoeff = permute(unmixC, [2 1 3 4]);
        kernelIm = permute(kernelIm, [2 1 3 4]);
        
        if ( f == 1 )
            gFactor = permute(gMex, [2 1]) / R;
        end
    
        disp(['kernel computation time is ' num2str(toc)]);
    end 
    
    if ( ~useACSLines )
        im = HTGRAPPA_SrcDstChannels_ImagedomainReconWithUnmixCoeff(reduced_k_data, unmixCoeff, zeroFilledSize);
        imfft = fft2DForVolume(im);
        imfft = performRawDataFilter(imfft, rawFilterFE, rawFilterPE);
        fullkspace(:,:,f) = imfft;
        unwrappedIm(:,:,f) = im;
    else
        option.acqPELines = sampledLineLoc(:,f);
        
%         fkspace = GRAPPA_SrcDstChannels_Recon_2D(reduced_k_data, Coef, option);
        
        aliasedIm = ifft2c(reduced_k_data);
        aliasedIm = repmat(aliasedIm, [1 1 1 dstCHA]);
        unWarppedIm = sum(aliasedIm.*kernelIm, 3);
        unWarppedIm = squeeze(unWarppedIm);

        fkspace = fft2c(unWarppedIm);
        
        fkspace(:,sampling_location_ref(:),:) = oriRef(:,sampling_location_ref(:),:,f);
        im = ifft2c(fkspace);
        if ( ~isempty(zeroFilledSize) )
            im = Zero_Padding_Resize_NoFiltering(im, zeroFilledSize(1), zeroFilledSize(2));
        end
        fullkspace(:,:,:,f) = fft2c(im);
        unwrappedIm(:,:,:,f) = im;
    end

    sensitivityMap(:,:,:,f) = senMap(:,:,:,f);
    disp(['recon time is ' num2str(toc)]);
end

disp(['Total recon time is ' num2str(toc(tstart))]);
