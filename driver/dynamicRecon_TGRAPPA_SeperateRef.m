
function [unwrappedIm, fullkspace, sensitivityMap] = dynamicRecon_TGRAPPA_SeperateRef(underSampledKspace, Psi, ref, reductionFactor, options)
% This function performs the dynamic image reconstruction using TGRAPPA
% help info : please execute "help produceHelpMrRecon"

% ===============================================================
% prepare the MrRecon

s = size(underSampledKspace);
Nfe = s(1);
Npe = s(2);
numOfCoil = s(3);
numOfFrames = s(4);

Nfe_ref = size(ref,1);
Npe_ref = size(ref,2);
numOfCoil_ref = size(ref,3);
numOfFrames_ref = size(ref,4);

reconStrategy = lower(options.reconStrategy);
sensitivity = options.sensitivity;
if ( ~isempty(sensitivity) )    
    numOfSen = size(sensitivity, 4);
    if ( numOfSen < numOfFrames )
        error('At least one sensitivity must be available for every frame ... ');
    end
end

reference = lower(options.reference);
subtractMean = options.subtractMean;

% ----------------------------
% for GRAPPA and HTGRAPPA
linesUsedForGRAPPA = lower(options.linesUsedForGRAPPA);
updateRef = lower(options.updateRef);
minFEUsed = options.minFEUsed;
maxFEUsed = options.maxFEUsed;
blockCoord= options.blockCoord;
missingLines = options.missingLines;
firstDataLineInBlcok = options.firstDataLineInBlcok;
kernalIndFE = options.kernalIndFE;
useACSLines = options.useACSLines;
regularization = options.regularization;
thresReg = options.thresReg;
overDetermineRatio = options.overDetermineRatio;
kspaceCenterFE = options.kspaceCenterFE;
phaseEncodingUsed = options.phaseEncodingUsed;
voxelsize = options.voxelsize;
iters = options.iters;
sigma = options.sigma;

numOfBlocksForRef = options.numOfBlocksForRef;
timing = options.timing;
spatialSmoothingKernel = 32;

% ----------------------------
% detect the sampling lines
sampledLineLoc = detectSampledLinesDynamic(underSampledKspace);
sampledLineLoc_ref = detectSampledLinesDynamic(ref);
% ----------------------------
% detect the valide range for FE sampled data
sampledRangeFE = detectSampledRangeFE(underSampledKspace);
sampledRangeFE_ref = detectSampledRangeFE(ref);
% ----------------------------
% generate the filter if required
if ( strcmp(options.rawFilterFE, 'None') == 0 )
    rawFilterFE = Matlab_ComputeKSpaceFilter(Nfe, Nfe/2, sampledRangeFE(1), sampledRangeFE(2)-sampledRangeFE(1)+1, kspaceCenterFE, ...
                            options.rawFilterFE, options.rawFilterFEStrength, 0, 0);
    rawFilterFE = single(rawFilterFE);
else
    rawFilterFE = [];
end

if ( strcmp(options.rawFilterPE, 'None') == 0 )
    rawFilterPE = Matlab_ComputeKSpaceFilter(Npe, Npe/2, 0, Npe+1, Npe/2+1, ...
                            options.rawFilterPE, options.rawFilterPEStrength, 0, 0);
    rawFilterPE = single(rawFilterPE);
else
    rawFilterPE = [];
end

if ( strcmp(options.ascFilterFE, 'None') == 0 )
    ascFilterFE = Matlab_ComputeKSpaceFilter(Nfe_ref, Nfe_ref/2, sampledRangeFE_ref(1), sampledRangeFE_ref(2)-sampledRangeFE_ref(1)+1, kspaceCenterFE, ...
                            options.ascFilterFE, options.ascFilterFEStrength, 0, 0);
    ascFilterFE = single(ascFilterFE);
else
    ascFilterFE = [];
end

if ( strcmp(options.ascFilterPE, 'None') == 0 )
    ascFilterPE = Matlab_ComputeKSpaceFilter(Npe_ref, Npe_ref/2, 0, Npe_ref+1, Npe_ref/2+1, ...
                            options.ascFilterPE, options.ascFilterPEStrength, 0, 0);
    ascFilterPE = single(ascFilterPE);
else
    ascFilterPE = [];
end
% ----------------------------

% noise prewhitening
if ( options.preWhitening )
    for f=1:numOfFrames
        underSampledKspace(:,:,:,f) = PreWhittening( underSampledKspace(:,:,:,f), Psi );
    end
    
    for f=1:numOfFrames_ref
        ref(:,:,:,f) = PreWhittening( ref(:,:,:,f), Psi );
    end
end

% ===============================================================

% prepare the reference kspace

tstart = tic;

disp('TGRAPPA starting ... ');    

% allocate outputs
fullkspace = zeros(s, 'single');
unwrappedIm = zeros(s, 'single');
sensitivityMap = zeros(s, 'single');

% for f=1:numOfFrames_ref    
%     currRef = ref(:,:,:,f);    
%     currRef = performRawDataFilter(currRef, ascFilterFE, ascFilterPE);
%     complexImageSen = ifft2DForVolume(currRef);
%     complexImageSen = Zero_Padding_Resize_NoFiltering(complexImageSen, s(1), s(2));
%     [combinedImSen, senMap] = adaptiveCoilCombination2DForSensitivity_PKMethod(complexImageSen, spatialSmoothingKernel);
%     sensitivityMap(:,:,:,f) = senMap;
% end

for f=numOfFrames_ref+1:numOfFrames
    sensitivityMap(:,:,:,f) = sensitivityMap(:,:,:,numOfFrames_ref);
end

% subtract mean
if ( subtractMean )
    comp = computeTemporalMean(underSampledKspace, reductionFactor);
    compIm = ifft2DForVolume(comp);
end

sampledLineLoc_ref = repmat([1:reductionFactor:Npe_ref]', [1 numOfFrames]);

for f=1:numOfFrames

    if ( timing ) 
        disp(['Frame ' num2str(f)]);
        tic; 
    end

    % generate the header for current frame
    header = generateDynamicReconHeader(f, reductionFactor, blockCoord, Nfe_ref, Npe_ref, numOfCoil, minFEUsed, maxFEUsed, sampledLineLoc_ref);        
    % generate header for acs data
    headerAcs = generateDynamicReconHeaderACS(header, reconStrategy, linesUsedForGRAPPA, reductionFactor, Npe_ref);

    if ( numOfFrames > 1 )
        reduced_k_data = underSampledKspace(:,:,:,f);    
    else
        reduced_k_data = underSampledKspace(:,:,:);
    end
    
    if ( numOfFrames_ref > 1 )
        currRef = ref(:,:,:,f);
    else
        currRef = ref(:,:,:);
    end
    
    if ( ~isempty(ascFilterFE) )
        currRef = performRawDataFilter(currRef, ascFilterFE, ascFilterPE);
    end
    
    currRef = currRef(sampledRangeFE_ref(1):sampledRangeFE_ref(2), :, :);
    
    % grappa recon   
                        
%     option.KernelSize = [length(blockCoord) length(kernalIndFE)];
%     option.KernelPattern = blockCoord;
%     option.OutPattern = missingLines;
%     option.SamplingPattern = sampledLineLoc(:,f);
%     
%     Coef = GRAPPA_Kernel_2D(currRef, option, thresReg);       
%     fkspace = GRAPPA_Recon_2D(reduced_k_data, Coef, option);
    % plotKSpace(fkspace);
    
    % try my code
    option.KernelSize   = [length(kernalIndFE) length(blockCoord)];
    option.acqPELines   = 1:size(currRef,2);
    option.acqAcsLines  = 1:size(currRef,2);
    option.KernelPattern = blockCoord;
    option.OutPattern = missingLines;

    Coef = GRAPPA_SrcDstChannels_Kernel_2D(currRef, currRef, option, thresReg);
    
    option.acqPELines = sampledLineLoc(:,f);
    fkspace =  GRAPPA_SrcDstChannels_Recon_2D(reduced_k_data, Coef, option);
    % plotKSpace(fkspace);

    if ( options.useACSLines )
        fkspace(:,options.sampling_location_ref(:),:) = options.original_ref(:,options.sampling_location_ref(:),:, f);
    end
    
    if ( ~isempty(rawFilterFE) )
        fkspace = performRawDataFilter(fkspace, rawFilterFE, rawFilterPE);
    end
    
    if ( numOfFrames > 1 )
        fullkspace(:,:,:,f) = fkspace;
        unwrappedIm(:,:,:,f) = ifft2DForVolume(fkspace);
    else
        fullkspace = fkspace;
        unwrappedIm = ifft2DForVolume(fkspace);        
    end
    
    if ( timing ) 
        disp(['recon time is ' num2str(toc)]);
    end
end 

disp(['Total recon time is ' num2str(toc(tstart))]);
