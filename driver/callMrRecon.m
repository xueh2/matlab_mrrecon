
function [unwarppedImSenseResized, unwarppedImGrappaResized, unwarppedImLinearResized, unwarppedImVICPAASResized, fullKSpaceGrappa, fullKSpaceLinear, fullKSpaceVICPAAS, voxelsize, gFactorGrappa, gFactorSense] = callMrRecon(kspace, Noise, rx_dwelltime_noise, ref, phsCorr, reflect, reflectRef, reflectPhsCorr, protocol, headers, option, VBorVD)

% script for perform all dimensional recon

% <!-- EPI Parameters -->

% <p><s>MEAS.sPat.lAccelFactPE</s>        <d>gadgetron.epiparams.acceleration.value</d></p>
% <p><s>MEAS.sPat.lRefLinesPE</s>         <d>gadgetron.epiparams.accreflines.value</d></p>
% <p><s>YAPS.alRegridMode</s>             <d>gadgetron.epiparams.rampsampling.value</d></p>
% <p><s>YAPS.alRegridRampupTime</s>       <d>gadgetron.epiparams.rampuptime.value</d></p>
% <p><s>YAPS.alRegridFlattopTime</s>      <d>gadgetron.epiparams.flattoptime.value</d></p>
% <p><s>YAPS.alRegridRampdownTime</s>     <d>gadgetron.epiparams.rampdowntime.value</d></p>
% <p><s>YAPS.alRegridDelaySamplesTime</s> <d>gadgetron.epiparams.acqdelay.value</d></p>
% <p><s>YAPS.aflRegridADCDuration</s>     <d>gadgetron.epiparams.readouttime.value</d></p>
% <p><s>YAPS.alRegridDestSamples</s>      <d>gadgetron.epiparams.numsamples.value</d></p>
% <p><s>YAPS.lEchoSpacing</s>             <d>gadgetron.epiparams.echospacing.value</d></p>

% VBorVD : 1, VB; 2, VD

isVD = (VBorVD==2);

%% incoming kspace
s = size(kspace)
size(Noise)
sRef = size(ref)
sPhsCorr = size(phsCorr)

COL = s(1);
LIN = s(2);
CHA = s(3);

unwarppedImSenseResized = [];
unwarppedImGrappaResized = [];
unwarppedImLinearResized = [];
unwarppedImVICPAASResized = [];

fullKSpaceGrappa = [];
fullKSpaceLinear = [];
fullKSpaceVICPAAS = [];
gFactorGrappa = [];
gFactorSense = [];

%% find the pat mode
RefScanMode = protocol.sPat.ucRefScanMode;

if ( ~isVD )
    ind = find(RefScanMode=='0');
    RefScanMode = RefScanMode(ind:end)
end

FOV_reduction_factor = protocol.sPat.lAccelFactPE
FOV_reduction_factor3D = protocol.sPat.lAccelFact3D

if ( FOV_reduction_factor>1 | FOV_reduction_factor3D>1 )
    numRefLines = protocol.sPat.lRefLinesPE
else
    numRefLines = 0;
end

if ( isVD )
    rx_dwelltime_data = protocol.sRXSPEC.alDwellTime_{1};
else
    rx_dwelltime_data = protocol.sRXSPEC.alDwellTime{1};
end
[feFOV, peFOV, sliceThickness] = findFOVFromConfig(headers.Config)

if ( isVD )
    sliceThickness = protocol.sSliceArray.asSlice_{1}.dThickness;
end

performSENSE = option.performSENSE;
performGRAPPA = option.performGRAPPA;
performVICPAAS = option.performVICPAAS;

nopat = 0;
ipat = 0;
spat = 0;
tpat = 0;
avepat = 0;

if ( isVD )
    if ( RefScanMode==1 )
        nopat = 1;
    end
    
    if ( RefScanMode==2 )
        ipat = 1;
    end

    if ( RefScanMode==4 | RefScanMode==8 )
        spat = 1;
    end

    if ( RefScanMode==32 | RefScanMode==64 )
        tpat = 1;
    end

    if ( RefScanMode==16 )
        avepat = 1;
    end    
else
    
%   enum PATRefScanMode
%   {
%     PAT_REF_SCAN_UNDEFINED      = 0x01, // e.g. if no PAT is selected
%     PAT_REF_SCAN_INPLACE        = 0x02, // sequence supplies inplace reference lines
%     PAT_REF_SCAN_EXTRA          = 0x04, // sequence supplies extra reference lines
%     PAT_REF_SCAN_PRESCAN        = 0x08, // sequence does not supply reference lines, the data must have been acquired with a previous measurement
%     PAT_REF_SCAN_INTRINSIC_AVE  = 0x10, // The sequence contains intrinsic ref.lines due to sharing e.g. in the averages dimension
%     PAT_REF_SCAN_INTRINSIC_REP  = 0x20, // The sequence contains intrinsic ref.lines due to sharing e.g. in the repetition or phases dimension (i.e., TSENSE)
%     PAT_REF_SCAN_INTRINSIC_PHS  = 0x40, // The sequence contains intrinsic ref.lines due to sharing e.g. in the repetition or phases dimension (i.e., TSENSE)
%     PAT_REF_SCAN_INPLACE_LET    = 0x80  // A single (L)ong (E)cho (T)rain acquires reference lines and imaging lines
%   };

    if ( ~isempty(strfind(RefScanMode, '0x1')) & isempty(strfind(RefScanMode, '0x10')) )
        nopat = 1;
    end
    
    if ( ~isempty(strfind(RefScanMode, '0x4')) & isempty(strfind(RefScanMode, '0x40')) )
        spat = 1;
    end

    if ( ~isempty(strfind(RefScanMode, '0x2')) & isempty(strfind(RefScanMode, '0x20')) )
        ipat = 1;
    end

    if ( ~isempty(strfind(RefScanMode, '0x20')) | ~isempty(strfind(RefScanMode, '0x40')) )
        tpat = 1;
    end

    if ( ~isempty(strfind(RefScanMode, '0x10')) )
        avepat = 1;
    end
end

if ( ipat )    
    if ( size(ref,1) ~= size(kspace, 1) )
        ipat = 0;
        spat = 1;
    end
end

nopat
ipat
spat
tpat
avepat

%% get the incoming ICE dimension
% [Col Line Cha Acq Slice Partition Echo Phase Rep Set Seg]
ACQ = 1;
SLC = 1;
PAR = 1;
ECO = 1;
PHS = 1;
REP = 1;
SET = 1;
SEG = 1;

maxDimKSpace = 3;

if ( numel(s) >= 4 ) 
    ACQ = s(4); 
end

if ( numel(s) >= 5 ) 
    SLC = s(5); 
end

if ( numel(s) >= 6 ) PAR = s(6); end
if ( numel(s) >= 7 ) ECO = s(7); end
if ( numel(s) >= 8 ) PHS = s(8); end
if ( numel(s) >= 9 ) REP = s(9); end
if ( numel(s) >= 10 ) SET = s(10); end
if ( numel(s) >= 11 ) SEG = s(11); end

if ( isempty(phsCorr) & (SEG>1) )
    kspace2 = sum(kspace, 11);
    kspace = kspace2;
    s = size(kspace)
end

maxDimKSpace = numel(s);

% ref
if ( ~isempty(ref) )
    refCOL = sRef(1);
    refLIN = sRef(2);
    refCHA = sRef(3);

    refACQ = 1;
    refSLC = 1;
    refPAR = 1;
    refECO = 1;
    refPHS = 1;
    refREP = 1;
    refSET = 1;
    refSEG = 1;

    maxDimRef = 3; 

    if ( numel(sRef) >= 4 ) refACQ = sRef(4); end
    if ( numel(sRef) >= 5 ) refSLC = sRef(5); end
    if ( numel(sRef) >= 6 ) refPAR = sRef(6); end
    if ( numel(sRef) >= 7 ) refECO = sRef(7); end
    if ( numel(sRef) >= 8 ) refPHS = sRef(8); end
    if ( numel(sRef) >= 9 ) refREP = sRef(9); end
    if ( numel(sRef) >= 10 ) refSET = sRef(10); end
    if ( numel(sRef) >= 11 ) refSEG = sRef(11); end

    if ( isempty(phsCorr) & (refSEG>1) )
        ref2 = sum(ref, 11);
        ref = ref2;
        sRef = size(ref)
    end

    maxDimRef = numel(sRef);
end

% ref
if ( ~isempty(phsCorr) )
    phsCorrCOL = sPhsCorr(1);
    phsCorrLIN = sPhsCorr(2);
    phsCorrCHA = sPhsCorr(3);

    phsCorrACQ = 1;
    phsCorrSLC = 1;
    phsCorrPAR = 1;
    phsCorrECO = 1;
    phsCorrPHS = 1;
    phsCorrREP = 1;
    phsCorrSET = 1;
    phsCorrSEG = 1;

    maxDimphsCorr = 3; 

    if ( numel(sPhsCorr) >= 4 ) phsCorrACQ = sPhsCorr(4); end
    if ( numel(sPhsCorr) >= 5 ) phsCorrSLC = sPhsCorr(5); end
    if ( numel(sPhsCorr) >= 6 ) phsCorrPAR = sPhsCorr(6); end
    if ( numel(sPhsCorr) >= 7 ) phsCorrECO = sPhsCorr(7); end
    if ( numel(sPhsCorr) >= 8 ) phsCorrPHS = sPhsCorr(8); end
    if ( numel(sPhsCorr) >= 9 ) phsCorrREP = sPhsCorr(9); end
    if ( numel(sPhsCorr) >= 10 ) phsCorrSET = sPhsCorr(10); end
    if ( numel(sPhsCorr) >= 11 ) phsCorrSEG = sPhsCorr(11); end

    maxDimPhsCorr = numel(sPhsCorr);
end

%% if it is with ramp sampling, correct the kspace


%% remove read-out oversampling
kspace = performDownSampleFE(kspace);
COL = size(kspace, 1);
size(kspace)

ref = performDownSampleFE(ref);
refCOL = size(ref, 1);
size(ref)

phsCorr = performDownSampleFE(phsCorr);

%% if need phase correction, do it
if ( ~isempty(phsCorr) )
    
    % flip the phase corr data
    for seg=1:phsCorrSEG
        for set=1:phsCorrSET
            for rep=1:phsCorrREP
                for phs=1:phsCorrPHS
                    for eco=1:phsCorrECO
                        for par=1:phsCorrPAR
                            for slc=1:phsCorrSLC                        
                                for acq=1:phsCorrACQ
                                    for lin=1:phsCorrLIN
                                        needFlip = reflectPhsCorr(1, lin, 1, acq, slc, par, eco, phs, rep, set, seg);
                                        if ( needFlip )
                                            phsData = phsCorr(:, lin, :, acq, slc, par, eco, phs, rep, set, seg);
                                            phsData = flipdim(phsData, 1);
                                            phsCorr(:, lin, :, acq, slc, par, eco, phs, rep, set, seg) = phsData;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % flip the kspace
    for seg=1:SEG
        for set=1:SET
            for rep=1:REP
                for phs=1:PHS
                    for eco=1:ECO
                        for par=1:PAR
                            for slc=1:SLC                        
                                for acq=1:ACQ
                                    for lin=1:LIN
                                        needFlip = reflect(1, lin, 1, acq, slc, par, eco, phs, rep, set, seg);
                                        if ( needFlip )
                                            kspaceData = kspace(:, lin, :, acq, slc, par, eco, phs, rep, set, seg);
                                            kspaceData = flipdim(kspaceData, 1);
                                            kspace(:, lin, :, acq, slc, par, eco, phs, rep, set, seg) = kspaceData;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % flip the ref
    if ( ~isempty(ref) & ~isempty(reflectRef) )
        for seg=1:refSEG
            for set=1:refSET
                for rep=1:refREP
                    for phs=1:refPHS
                        for eco=1:refECO
                            for par=1:refPAR
                                for slc=1:refSLC                        
                                    for acq=1:refACQ
                                        for lin=1:refLIN
                                            needFlip = reflectRef(1, lin, 1, acq, slc, par, eco, phs, rep, set, seg);
                                            if ( needFlip )
                                                refData = ref(:, lin, :, acq, slc, par, eco, phs, rep, set, seg);
                                                refData = flipdim(refData, 1);
                                                ref(:, lin, :, acq, slc, par, eco, phs, rep, set, seg) = refData;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % perform phase correction
    
    if ( option.PerformPhaseCorrection )
        for set=1:SET
            for eco=1:ECO
                for par=1:PAR
                    for slc=1:SLC 
                        for rep=1:REP
                            for phs=1:PHS

                                reduced_k_data = kspace(:,:,:,:, slc, par, eco, phs, rep, set, :);                    
                                reduced_k_data = permute(reduced_k_data, [1 2 3 4 11 5 6 7 8 9 10]);
                                reduced_k_data = squeeze(reduced_k_data);                            

                                reduced_k_data = reshape(reduced_k_data, [COL LIN CHA 1 SEG 1 1]);

                                if ( phsCorrREP>1 )
                                    phaseDataSlice = phsCorr(:, phsCorrLIN/2+1, :, :, slc, par, eco, 1, rep, set, :);
                                elseif ( phsCorrPHS>1 )
                                    phaseDataSlice = phsCorr(:, phsCorrLIN/2+1, :, :, slc, par, eco, phs, 1, set, :);
                                else
                                    phaseDataSlice = phsCorr(:, phsCorrLIN/2+1, :, :, slc, par, eco, 1, 1, set, :);
                                end

                                phaseDataSlice = squeeze(phaseDataSlice);
%                                 phaseDataSlice(:,:,1,2) =repmat(sum( conj(phaseDataSlice(:,:,1,1)).*phaseDataSlice(:,:,1,2),2),[1 CHA]);
%                                 phaseDataSlice(:,:,2,2) =repmat(sum( conj(phaseDataSlice(:,:,1,1)).*phaseDataSlice(:,:,2,2),2),[1 CHA]);
%                                 phaseDataSlice(:,:,1,1) = ones(COL,CHA);

                                if ( isempty(reflectRef) )
                                    phaseDataSliceSum = sum(phaseDataSlice, 2);
                                    for c=1:CHA
                                        % phaseDataSlice(:,c,:,:) = phaseDataSlice(:,1,:,:);
                                        phaseDataSlice(:,c,:,:) = phaseDataSliceSum;
                                    end
                                end
                                
                                filterPhaseData = 1;
                                robustEstimation = 0;
                                reduced_k_data = performEPIPhaseCorrection(reduced_k_data, phaseDataSlice, filterPhaseData, robustEstimation, option.PhaseCorrMethod);
                                % reduced_k_data = performEPIPhaseCorrection(reduced_k_data, phaseDataSlice, filterPhaseData, robustEstimation, 'LinearPhase');

                                reduced_k_data = squeeze(reduced_k_data);
                                reduced_k_data = permute(reduced_k_data, [1 2 3 5 4]);
                                kspace(:,:,:,:, slc, par, eco, phs, rep, set, :) = reduced_k_data;

                                if ( ~isempty(ref) & spat & ~isempty(reflectRef) )
                                    if ( (rep<=refREP) & (phs<=refPHS) )
                                        refSlice = ref(:,:,:,:, slc, par, eco, phs, rep, set, :);
                                        refSlice = permute(refSlice, [1 2 3 4 11 5 6 7 8 9 10]);
                                        refSlice = squeeze(refSlice);
                                        refSlice = reshape(refSlice, [refCOL refLIN refCHA 1 refSEG 1 1]);
                                        refSlice = performEPIPhaseCorrection(refSlice, phaseDataSlice, filterPhaseData, robustEstimation, option.PhaseCorrMethod);
                                        refSlice = squeeze(refSlice);
                                        refSlice = permute(refSlice, [1 2 3 5 4]);
                                        ref(:,:,:,:, slc, par, eco, phs, rep, set, :) = refSlice;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % sum over SEG
    if ( SEG > 1 )
        kspace2 = sum(kspace, 11);
        kspace = kspace2;
        s = size(kspace)
    end

    if ( ~isempty(ref) & refSEG>1 )
        ref2 = sum(ref, 11);
        ref = ref2;
        sRef = size(ref)
    end    
end

%% noise prewhitening
disp('performing noise prewhitening...') 

if ( ~isempty(Noise) )
    Noise = squeeze(Noise);
    Noise = permute(Noise, [2 1 3]);
    noisePrewhiteningMatrix = calculateNoisePrewhitener(Noise, rx_dwelltime_data, rx_dwelltime_noise);

    kspace2 = permute(kspace, [3 1 2 4:maxDimKSpace]);
    kspace2 = applyNoisePrewhitener(kspace2, noisePrewhiteningMatrix);
    kspace2 = permute(kspace2, [2 3 1 4:maxDimKSpace]);
    size(kspace2)
    kspace = kspace2;
    size(kspace)
    clear kspace2

    if ( ~isempty(ref) )
        ref2 = permute(ref, [3 1 2 4:maxDimRef]);
        ref2 = applyNoisePrewhitener(ref2, noisePrewhiteningMatrix);
        ref2 = permute(ref2, [2 3 1 4:maxDimRef]);
        size(ref2)
        ref = ref2;
        size(ref)
        clear ref2
    end
end

newNpe = round(peFOV/(feFOV/COL))
if ( LIN > newNpe )
    newNpe = LIN;
end

newNpar = protocol.sKSpace.lPartitions;

voxelsize = [feFOV/COL peFOV/newNpe sliceThickness/newNpar]

% plotKSpaceArray(kspace);

%% check whether it is the random sampling, just for the CS support
isRandomSampling = 0;
if ( nopat )
    
    kspace2D = kspace(:,:,1,1,1,1,1,1,1,1);
    sampledLineLoc = detectSampledLines(kspace2D);
    
    samplingStepSize = sampledLineLoc(2:end) - sampledLineLoc(1:end-1);
    
    uniqueSampledStepSize = unique(samplingStepSize);
    if ( numel(uniqueSampledStepSize) > 1 )
        isRandomSampling = 1;
        disp('Incoming scan is randomly sampled ... ');        
        tpat = 1;
    end
end

%% prepare the kspace and ref
% kspace will be [COL LIN CHA PAR ]

if ( tpat )
    
    tpat_useREP = 0;
    tpat_usePHS = 0;
    
    if ( isVD )
        % REP used as tpat dim
        if ( RefScanMode == 32 )
            tpat_useREP = 1;
        end

        % PHS used as tpat dim
        if ( RefScanMode == 64 )
            tpat_usePHS = 1;
        end        
    else
        % REP used as tpat dim
        if ( strfind(RefScanMode, '0x20') )
            tpat_useREP = 1;
        end

        % PHS used as tpat dim
        if ( strfind(RefScanMode, '0x40') )
            tpat_usePHS = 1;
        end
    end
    
    if ( isRandomSampling )
        % incoming scan may be a tpat scan
        if ( PHS>1 & REP==1 & PAR==1 )
            tpat_usePHS = 1;
            tpat_useREP = 0;
        end

        if ( REP>1 & PHS==1 & PAR==1 )
            tpat_useREP = 1;
            tpat_usePHS = 0;
        end
    end
    
    option_ori = option;
    option.KernelSizeSPIRiT = option.KernelSizeSPIRiT_TPAT;
    option.OutKernelSPIRiT = option.OutKernelSPIRiT_TPAT;
    option.thresRegSPIRiT = option.thresRegSPIRiT_TPAT;
end

if ( spat )
    if ( FOV_reduction_factor3D == 1 )
        if ( refPAR == 1 )
            loc = detectSampledLines(ref(:,:,1,1,1,1,1,1,1,1))
        else
            loc = detectSampledLines(ref(:,:,1,1,1,round(refPAR/2),1,1,1,1))
        end

        stepSize = loc(2:end) - loc(1:end-1);

        if ( ~tpat )
            if ( ~isempty(find(stepSize>1)) )
                spat = 0;
                ipat = 1;
            else    
                ref= ref(:,loc(1):loc(end), :,:,:,:,:,:,:,:);
            end
        end
    else
        sumRef = sum(ref, 6);        
        locLIN = detectSampledLines(squeeze(sumRef))
        
        sumPAR = squeeze(sum(ref, 2));
        sumPAR = permute(sumPAR, [1 3 2 4]);
        locPAR = detectSampledLines(squeeze(sumPAR))
        
        ref= ref(:,locLIN(1):locLIN(end), :,:,:,locPAR(1):locPAR(end),:,:,:,:);
    end
end

if ( ipat )          
    kspace_ipat_withRef = kspace;

    ref_ipat = ref;
    
    if ( FOV_reduction_factor3D == 1 )
        loc = detectSampledLines(ref_ipat(:,:,1,1,1,1,1,1,1,1))
        diff = loc([2:end 1]) - loc;
        ind = find(diff==1);
        startLine = loc(ind(1));
        endLine = loc(ind(end))+1;
        ref_ipat= ref_ipat(:,startLine:endLine, :,:,:,:,:,:,:,:);
        
        kspace_ipat_withRef(:,startLine:endLine,:,:,:,:,:,:,:,:) = ref_ipat;
    else
        sumRef = sum(ref_ipat, 6);        
        locLIN = detectSampledLines(squeeze(sumRef))
        
        sumPAR = squeeze(sum(ref_ipat, 2));
        sumPAR = permute(sumPAR, [1 3 2 4]);
        locPAR = detectSampledLines(squeeze(sumPAR))
        
        ref_ipat= ref_ipat(:,locLIN(1):locLIN(end), :,:,:,locPAR(1):locPAR(end),:,:,:,:);
        
        kspace_ipat_withRef(:,locLIN(1):locLIN(end),:,:,:,locPAR(1):locPAR(end),:,:,:,:) = ref_ipat;
    end
    % plotKSpaceSamplingPattern(ref_ipat(:,:,1,:,:,:,:,:,:,:)); 
end

if ( avepat )
    % if  ACQ > 1 and ref is empty, then sum up ACQ as reference
    if ( ACQ>1 & isempty(ref) )
        ref = sum(kspace, 4);
        
        % only the central kspace shall be used for acs
        startLIN = LIN/2 - numRefLines/2 + 1;
        endLIN = startLIN + numRefLines;        
        ref = ref(:,startLIN:endLIN, :,:,:);
    end
end

%% handle the partial fourier and asymmetric echo properly
SamplesInScan = option.SamplesInScan;
KSpaceCentreColumn = option.KSpaceCentreColumn;
MaxKSpaceLineNo = option.MaxKSpaceLineNo;
KSpaceCentreLineNo = option.KSpaceCentreLineNo;
KSpaceCentrePartitionNo = option.KSpaceCentrePartitionNo;

% LIN
while ( mod(MaxKSpaceLineNo, FOV_reduction_factor)~=0 )
    MaxKSpaceLineNo = MaxKSpaceLineNo + 1;
end

rangeUsedLIN = [];
aZ = addPrePostZeros(KSpaceCentreLineNo, LIN);

if ( aZ==1 )
    rangeUsedLIN = [LIN-MaxKSpaceLineNo+1 LIN];
end

if ( aZ==2 )
    rangeUsedLIN = [1 MaxKSpaceLineNo];
end

option.rangeUsedLIN = rangeUsedLIN;

if ( tpat )
    if ( size(kspace, 2) ~= MaxKSpaceLineNo )
        stpat = s;
        stpat(1) = size(kspace,1);
        stpat(2) = MaxKSpaceLineNo;
        kspace2 = zeros(stpat);
        kspace2(:,1:size(kspace,2),:,:,:,:,:,:,:,:,:) = kspace;
        kspace = kspace2;
        clear kspace2;
    end
end

% FE
rangeUsedFE = [];
if ( tpat | ipat | nopat )
    aZ = addPrePostZeros(KSpaceCentreColumn, SamplesInScan);
    
    if ( aZ==1 )
        rangeUsedFE = [COL-SamplesInScan/2+1 COL];
    end
    
    if ( aZ==2 )
        rangeUsedFE = [1 SamplesInScan/2];
    end
end
option.rangeUsedFE = rangeUsedFE;
    
%% handle the partial fourier along the PAR
MaxKSpacePARNo = PAR;
while ( mod(MaxKSpacePARNo, FOV_reduction_factor3D)~=0 )
    MaxKSpacePARNo = MaxKSpacePARNo + 1;
end

startPAR = 1;
endPAR = PAR;
if ( FOV_reduction_factor3D > 1 )
    
    aZ = addPrePostZeros(KSpaceCentrePartitionNo, MaxKSpacePARNo);
    
    if ( aZ~=0 & ~strcmp(protocol.sKSpace.ucSlicePartialFourier, '0x10') )
        
        if ( aZ==1 )
            startPAR = newNpar - MaxKSpacePARNo + 1;
            endPAR = startPAR+PAR-1;
        end
        
        if ( aZ==2 )
            startPAR = 1;
            endPAR = PAR;
        end
               
        % extend kspace
        newKSpace = zeros(COL, LIN, CHA, ACQ, SLC, newNpar, ECO, PHS, REP, SET);
        newKSpace(:,:,:,:,:,startPAR:endPAR,:,:,:,:) = kspace(:,:,:,:,:,1:PAR,:,:,:,:);
        kspace = newKSpace;

        if ( ipat )
            newKSpace = zeros(COL, LIN, CHA, ACQ, SLC, newNpar, ECO, PHS, REP, SET);
            if ( endPAR < size(kspace_ipat_withRef, 6) )
                endPAR = size(kspace_ipat_withRef, 6);
            end
            newKSpace(:,:,:,:,:,startPAR:endPAR,:,:,:,:) = kspace_ipat_withRef;
            kspace_ipat_withRef = newKSpace;

            if ( aZ==1 )
                newRefPAR = refPAR+newNpar-MaxKSpacePARNo;
                newRef = zeros(refCOL, refLIN, refCHA, refACQ, refSLC, newRefPAR, refECO, refPHS, refREP, refSET);
                newRef(:,:,:,:,:,end-refPAR+1:end,:,:,:,:) = ref;
                ref = newRef;
            end
        end

        clear newKSpace newRef
    end
end
option.rangeUsedPAR = [startPAR endPAR];
        
PAR = size(kspace, 6);

%% compute the kspace filter

% image filter, always symmetric
filterCOL = generateKSpaceFilter(option.rawFilterFE, option.rawFilterFEStrength, COL, [1 COL], floor(COL/2));
filterLIN = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, LIN, [1 LIN], floor(LIN/2));
filterPAR = generateKSpaceFilter(option.rawFilterPAR, option.rawFilterPARStrength, PAR, [1 PAR], floor(PAR/2));

% ref coil map filter, symmetric along the sampled ref region
if ( ~isempty(ref) )
    if ( ~spat & ~isempty(option.rangeUsedFE) )
        range = findSymmetricSampledRegion(option.rangeUsedFE, floor(refCOL/2)+1)
    else
        range = [1 refCOL];
    end
    option.filterRefCOL = generateKSpaceFilter(option.acsFilterFE, option.acsFilterFEStrength, refCOL, range, floor(refCOL/2));

    rangePE = detectSampledRangePE(ref);
    range = findSymmetricSampledRegion(rangePE, floor(size(ref,2)/2)+1)
    option.filterRefLIN = generateKSpaceFilter(option.acsFilterPE, option.acsFilterPEStrength, size(ref,2), range, floor(size(ref,2)/2));

    if ( refPAR > 1 )
        refPAR = sum(ref, 3);
        refPAR = sum(refPAR, 4);
        refPAR = sum(refPAR, 5);
        refPAR = sum(refPAR, 1);
        refPAR = refPAR(1,:,1,1,1,:,1,1,1,1,1,1,1);
        refPAR = squeeze(refPAR);
        rangePAR = detectSampledRangePE(refPAR);

        range = findSymmetricSampledRegion(rangePAR, floor(PAR/2)+1)
        option.filterRefPAR = generateKSpaceFilter(option.acsFilterPE, option.acsFilterPEStrength, PAR, rangePAR, floor(PAR/2));
    else
        option.filterRefPAR = [];
    end
else
    option.filterRefCOL = filterCOL;
    option.filterRefLIN = filterLIN;
    option.filterRefPAR = filterPAR;
end

%% 2D or 2D + T recon
if ( nopat & ~isRandomSampling )
    % [Col Line Cha Acq Slice Partition Echo Phase Rep Set Seg]
    
    fullKSpaceGrappa = kspace;
    fullKSpaceLinear = kspace;
    fullKSpaceVICPAAS = kspace;
    
    unwarppedImGrappa = zeros(COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET);
        
    if ( PAR > 1 )
        kspace = ifftc(kspace,6);
    end
    
    for set=1:SET
        for rep=1:REP
            for phs=1:PHS
                for eco=1:ECO
                    for par=1:PAR
                        for slc=1:SLC                        
                            for acq=1:ACQ                                
                                kspaceCurr = kspace(:,:,:,acq, slc, par, eco, phs, rep, set);
                                senMap = CoilSensitivityMapEstimation(ifft2c(kspaceCurr), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
                                unwarppedImGrappa(:,:,acq, slc, par, eco, phs, rep, set) = SensitivityCoilCombination(ifft2c(kspaceCurr), senMap);
                            end
                        end
                    end
                end
            end
        end
    end
       
    if ( FOV_reduction_factor3D == 1 )
        unwarppedImGrappa = performKSpaceFilter2D(fft2c(unwarppedImGrappa), filterCOL, filterLIN);
        unwarppedImGrappa = ifft2c(unwarppedImGrappa);
        unwarppedImLinearResized = Zero_Padding_Resize_NoFiltering(unwarppedImGrappa, COL, newNpe);
    else
        unwarppedImGrappa = performKSpaceFilter3D(fft3c_ImageDimensions(unwarppedImGrappa), filterCOL, filterLIN, filterPAR);
        unwarppedImGrappa = ifft3c_ImageDimensions(unwarppedImGrappa);
        unwarppedImLinearResized = Zero_Padding_Resize_NoFiltering3D(unwarppedImGrappa, COL, newNpe, newNpar);
    end    
    
    unwarppedImSenseResized = unwarppedImLinearResized;
    unwarppedImGrappaResized = unwarppedImLinearResized;
    unwarppedImVICPAASResized = unwarppedImLinearResized;
else

if ( FOV_reduction_factor3D == 1 ) 
    
    if ( PAR > 1 )
        kspace = ifftc(kspace,6);
        ref = ifftc(ref,6);
        
        if ( ipat )
            kspace_ipat_withRef = ifftc(kspace_ipat_withRef,6);
            ref_ipat = ifftc(ref_ipat,6);
        end
        
        option.temporalScalingFactorSPIRiT = 2;
        option.dstChaThres  = 0.002;
        option.maxIterSPIRiT = 70;
        option.stopThresSPIRiT = 1e-5;
    end
    
    %% linear
    if ( performSENSE & ~isRandomSampling )
        
        % option.zeroFilledSize = [COL newNpe];
        
        if ( mod(LIN, FOV_reduction_factor) ~= 0 )
            LIN = LIN + FOV_reduction_factor - mod(LIN, FOV_reduction_factor);
            
            filterLIN = generateKSpaceFilter(option.rawFilterPE, option.rawFilterPEStrength, LIN, [0 LIN-1], floor(LIN/2));
            
            if ( ~empty(option.rangeUsedLIN) )
                range = findSymmetricSampledRegion(option.rangeUsedLIN, floor(LIN/2)+1)
            else
                range = [1 LIN];
            end
            option.filterRefLIN = generateKSpaceFilter(option.acsFilterPE, option.acsFilterPEStrength, LIN, range, floor(LIN/2));
        end
        
        unwarppedImSense = zeros(COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET);
        gFactorSense = zeros(COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET);       
        
        for set=1:SET
            for eco=1:ECO
                for slc=1:SLC
                    for par=1:PAR
                        for acq=1:ACQ

                            disp(['acq - par - slc - eco - set : ' num2str(acq) ' ' num2str(par) ' ' num2str(slc) ' ' num2str(eco) ' ' num2str(set)]);

                            if ( tpat )

                                if ( tpat_useREP )
                                    kspaceForRecon = kspace(:,:,:,acq, slc, par, eco, 1, :, set);
                                end

                                if ( tpat_usePHS )
                                    kspaceForRecon = kspace(:,:,:,acq, slc, par, eco, :, 1, set);
                                end
                               
                                kspaceForRecon = squeeze(kspaceForRecon);
                                [im, fullkspace, sensitivityMap, gFactor, im_oneSenMap, im_twoSenMap] = TSENSE_AverageAll_SNRUnit(kspaceForRecon, FOV_reduction_factor, option);
                                   
                                meanKSpace = mean(fullkspace, 4);
                                senMap = CoilSensitivityMapEstimation(ifft2c(meanKSpace), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.spatialSmoothingKernel);
                                im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);
                                            
                                if ( tpat_useREP )
                                    for rep=1:REP
                                        unwarppedImSense(:,:,acq, slc, par, eco, 1, rep, set) = im(:,:,rep);
                                        gFactorSense(:,:,acq, slc, par, eco, 1, rep, set) = gFactor;
                                    end
                                end

                                if ( tpat_usePHS )
                                    for phs=1:PHS
                                        unwarppedImSense(:,:,acq, slc, par, eco, phs, 1, set) = im(:,:,phs);
                                        gFactorSense(:,:,acq, slc, par, eco, phs, 1, set) = gFactor;
                                    end
                                end
                            end

                            if ( spat | avepat )                        
                                for rep=1:REP
                                    for phs=1:PHS
                                        kspaceForRecon = kspace(:,:,:,acq, slc, par, eco, phs, rep, set);
                                        if (avepat) 
                                            refForRecon = ref(:,:,:,1, slc, par, eco, phs, rep, set);
                                        else
                                            if ( refREP < rep )
                                                refForRecon = ref(:,:,:,acq, slc, par, eco, phs, refREP, set);
                                            else
                                                refForRecon = ref(:,:,:,acq, slc, par, eco, phs, rep, set);
                                            end
                                        end
                                        
                                        if ( refREP > 1 )                                                                                    
                                            [im, fullkspace, sensitivityMap, gFactor, im_oneSenMap, im_twoSenMap] = TSENSE_SeperateRef_SNRUnit(kspaceForRecon, refForRecon, FOV_reduction_factor, option);

                                            meanKSpace = mean(fullkspace, 4);
                                            senMap = CoilSensitivityMapEstimation(ifft2c(meanKSpace), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
                                            im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);
                                            
                                            unwarppedImSense(:,:,acq, slc, par, eco, phs, rep, set) = im;
                                            gFactorSense(:,:,acq, slc, par, eco, phs, rep, set) = gFactor;
                                        else
                                            kspaceForRecon = kspace(:,:,:,acq, slc, par, eco, phs, :, set);
                                            kspaceForRecon = squeeze(kspaceForRecon);
                                                                                       
                                            [im, fullkspace, sensitivityMap, gFactor, im_oneSenMap, im_twoSenMap] = TSENSE_SeperateRef_SNRUnit(kspaceForRecon, refForRecon, FOV_reduction_factor, option);
                                                 
                                            meanKSpace = mean(fullkspace, 4);
                                            senMap = CoilSensitivityMapEstimation(ifft2c(meanKSpace), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.spatialSmoothingKernel);
                                            im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);
                                            
                                            unwarppedImSense(:,:,acq, slc, par, eco, phs, 1:size(im,3), set) = im;
                                            gFactorSense(:,:,acq, slc, par, eco, phs, rep, set) = gFactor;
                                        end
                                    end
                                    
                                    if ( refREP == 1 )
                                        break;
                                    end
                                end
                            end

                            if ( ipat )
                                disp('ipat ... ');
                                for rep=1:REP
                                    for phs=1:PHS
                                        kspaceForRecon = kspace(:,:,:,acq, slc, par, eco, phs, rep, set);
                                        refForRecon = ref(:,:,:,acq, slc, par, eco, phs, rep, set);
                                        
                                        [im, fullkspace, sensitivityMap, gFactor, im_oneSenMap, im_twoSenMap, E_0, V] = TSENSE_InplaceRef_SNRUnit(kspaceForRecon, refForRecon, FOV_reduction_factor, option);
                                        
                                        unwarppedImSense(:,:,acq, slc, par, eco, phs, rep, set) = im_twoSenMap;
                                        gFactorSense(:,:,acq, slc, par, eco, phs, rep, set) = gFactor;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        unwarppedImSense = performKSpaceFilter2D(fft2c(unwarppedImSense), filterCOL, filterLIN);
        unwarppedImSense = ifft2c(unwarppedImSense);
        
        unwarppedImSenseResized = Zero_Padding_Resize_NoFiltering(unwarppedImSense, COL, newNpe);

        clear unwarppedImSense
    end
    
    if ( performGRAPPA & ~isRandomSampling )
        
        LIN = size(kspace, 2);
        
        fullKSpaceGrappa = kspace;              
        unwarppedImGrappa = zeros(COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET);
        gFactorGrappa = zeros(COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET);       
        
        for set=1:SET
            for eco=1:ECO
                for slc=1:SLC
                    for par=1:PAR
                        for acq=1:ACQ

                            disp(['acq - par - slc - eco - set : ' num2str(acq) ' ' num2str(par) ' ' num2str(slc) ' ' num2str(eco) ' ' num2str(set)]);

                            if ( tpat )

                                if ( tpat_useREP )
                                    kspaceForRecon = kspace(:,:,:,acq, slc, par, eco, 1, :, set);
                                end

                                if ( tpat_usePHS )
                                    kspaceForRecon = kspace(:,:,:,acq, slc, par, eco, :, 1, set);
                                end
                               
                                kspaceForRecon = squeeze(kspaceForRecon);
                                [im, fullkspace, sensitivityMap, gFactor] = TGRAPPA_AverageAll_SrcDstChannels_SNRUnit(kspaceForRecon, FOV_reduction_factor, option);
                                
                                meanKSpace = mean(fullkspace, 4);
                                senMap = CoilSensitivityMapEstimation(ifft2c(meanKSpace), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
                                im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);

                                im = PerformPartialFourierHandling2D(im, option, 1);
                                
                                dstCha = size(fullkspace, 3);
                                
                                if ( tpat_useREP )
                                    for rep=1:REP
                                        fullKSpaceGrappa(:,:,1:dstCha,acq, slc, par, eco, 1, rep, set) = fullkspace(:,:,:,rep);
                                        unwarppedImGrappa(:,:,acq, slc, par, eco, 1, rep, set) = im(:,:,rep);
                                        gFactorGrappa(:,:,acq, slc, par, eco, 1, rep, set) = gFactor;
                                    end
                                end

                                if ( tpat_usePHS )
                                    for phs=1:PHS
                                        fullKSpaceGrappa(:,:,1:dstCha,acq, slc, par, eco, phs, 1, set) = fullkspace(:,:,:,phs);
                                        unwarppedImGrappa(:,:,acq, slc, par, eco, phs, 1, set) = im(:,:,phs);
                                        gFactorGrappa(:,:,acq, slc, par, eco, phs, 1, set) = gFactor;
                                    end
                                end
                            end

                            if ( spat | avepat )                        
                                for rep=1:REP
                                    for phs=1:PHS
                                        kspaceForRecon = kspace(:,:,:,acq, slc, par, eco, phs, rep, set);
                                        
                                        setRef = set;
                                        if ( refSET < set )
                                            setRef = refSET;
                                        end
                                        
                                        if (avepat) 
                                            refForRecon = ref(:,:,:,1, slc, par, eco, phs, rep, setRef);
                                        else
                                            if ( refREP < rep )
                                                refForRecon = ref(:,:,:,acq, slc, par, eco, phs, refREP, setRef);
                                            else
                                                refForRecon = ref(:,:,:,acq, slc, par, eco, phs, rep, setRef);
                                            end
                                        end
                                        
                                        if ( norm(refForRecon(:)) == 0 )
                                            refForRecon = ref(:,:,:,acq, slc, par, eco, phs, rep, 1);
                                        end
                                        
                                        if ( refREP > 1 )                                        
                                            [im, fullkspace, sensitivityMap, gFactor] = TGRAPPA_SeperateRef_SrcDstChannels_SNRUnit(kspaceForRecon, refForRecon, FOV_reduction_factor, option);
                                            dstCha = size(fullkspace, 3);
                                            fullKSpaceGrappa(:,:,1:dstCha,acq, slc, par, eco, phs, rep, set) = fullkspace;

%                                             senMap = CoilSensitivityMapEstimation(ifft2c(fullkspace), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
%                                             im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);
                                            unwarppedImGrappa(:,:,acq, slc, par, eco, phs, rep, set) = im;
                                            gFactorGrappa(:,:,acq, slc, par, eco, phs, rep, set) = gFactor;
                                        else
                                            kspaceForRecon = kspace(:,:,:,acq, slc, par, eco, phs, :, set);
                                            kspaceForRecon = squeeze(kspaceForRecon);
                                            
                                            [im, fullkspace, sensitivityMap, gFactor] = TGRAPPA_SeperateRef_SrcDstChannels_SNRUnit(kspaceForRecon, refForRecon, FOV_reduction_factor, option);
                                                                                        
                                            dstCha = size(fullkspace, 3);
                                            fullKSpaceGrappa(:,:,1:dstCha,acq, slc, par, eco, phs, 1:size(fullkspace,4), set) = fullkspace;

%                                             senMap = CoilSensitivityMapEstimation(ifft2c(mean(fullkspace, 4)), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
%                                             im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);
                                            unwarppedImGrappa(:,:,acq, slc, par, eco, phs, 1:size(im,3), set) = im;
                                            gFactorGrappa(:,:,acq, slc, par, eco, phs, rep, set) = gFactor;
                                        end
                                    end
                                    
                                    if ( refREP == 1 )
                                        break;
                                    end
                                end
                            end

                            if ( ipat )
                                disp('ipat ... ');
                                for rep=1:REP
                                    for phs=1:PHS
                                        kspaceForRecon = kspace(:,:,:,acq, slc, par, eco, phs, rep, set);
                                        refForRecon = ref(:,:,:,acq, slc, par, eco, phs, rep, set);

                                        [im, fullkspace, sensitivityMap, gFactor, E_0, V] = TGRAPPA_InplaceRef_SrcDstChannels_SNRUnit(kspaceForRecon, refForRecon, FOV_reduction_factor, option);

                                        fullkspace = PerformPartialFourierHandling2D(ifft2c(fullkspace), option, 1); fullkspace = fft2c(fullkspace);
                                        
                                        dstCha = size(fullkspace, 3);
                                        fullKSpaceGrappa(:,:,1:dstCha,acq, slc, par, eco, phs, rep, set) = fullkspace;
                                        
                                        senMap = CoilSensitivityMapEstimation(ifft2c(fullkspace), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
                                        im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);                                                                               
                                        
                                        unwarppedImGrappa(:,:,acq, slc, par, eco, phs, rep, set) = im;
                                        gFactorGrappa(:,:,acq, slc, par, eco, phs, rep, set) = gFactor;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        unwarppedImGrappa = performKSpaceFilter2D(fft2c(unwarppedImGrappa), filterCOL, filterLIN);
        unwarppedImGrappa = ifft2c(unwarppedImGrappa);
        
        unwarppedImGrappaResized = Zero_Padding_Resize_NoFiltering(unwarppedImGrappa, COL, newNpe);

        clear unwarppedImGrappa
    end
    
    %% non-linear
    if ( performVICPAAS )
        
        LIN = size(kspace, 2);
        fullKSpaceLinear = kspace;
        unwarppedImLinear = zeros(COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET);

        fullKSpaceVICPAAS = kspace;
        unwarppedImVICPAAS = unwarppedImLinear;

        for set=1:SET
            for eco=1:ECO
                for acq=1:ACQ

                    if ( tpat )
                        for slc=1:SLC
                            if ( tpat_useREP )
                                kspaceForRecon = kspace(:,:,:,acq, slc, 1, eco, 1, :, set);
                                ref = computeMeanKSpace(squeeze(kspaceForRecon));
                            end

                            if ( tpat_usePHS )
                                kspaceForRecon = kspace(:,:,:,acq, slc, 1, eco, :, 1, set);
                                ref = computeMeanKSpace(squeeze(kspaceForRecon));
                            end

                            kspaceForRecon = squeeze(kspaceForRecon);                        

                            if ( option.tspirit & performGRAPPA & ~isRandomSampling )
                                if ( tpat_useREP )
                                    ref = fullKSpaceGrappa(:,:,:,acq, slc, 1, eco, 1, :, set);
                                end

                                if ( tpat_usePHS )
                                    ref = fullKSpaceGrappa(:,:,:,acq, slc, 1, eco, :, 1, set);
                                end
                            end

                            ref = squeeze(ref);
                            [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS(kspaceForRecon, ref, option);

                            dstCHA = size(kspaceLinear, 3);

                            if ( tpat_useREP )
                                for rep=1:REP
                                    fullKSpaceLinear(:,:,1:dstCHA,acq, slc, 1, eco, 1, rep, set) = kspaceLinear(:,:,:,rep);
                                    unwarppedImLinear(:,:,acq, slc, 1, eco, 1, rep, set) = unwrappedImVICPAASLinear(:,:,rep);
                                    fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, 1, eco, 1, rep, set) = kspaceVICPAAS(:,:,:,rep);
                                    unwarppedImVICPAAS(:,:,acq, slc, 1, eco, 1, rep, set) = unwrappedImVICPAAS(:,:,rep);
                                end
                            end

                            if ( tpat_usePHS )
                                for phs=1:PHS
                                    fullKSpaceLinear(:,:,1:dstCHA,acq, slc, 1, eco, phs, 1, set) = kspaceLinear(:,:,:,phs);
                                    unwarppedImLinear(:,:,acq, slc, 1, eco, phs, 1, set) = unwrappedImVICPAASLinear(:,:,phs);
                                    fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, 1, eco, phs, 1, set) = kspaceVICPAAS(:,:,:,phs);
                                    unwarppedImVICPAAS(:,:,acq, slc, 1, eco, phs, 1, set) = unwrappedImVICPAAS(:,:,phs);
                                end
                            end
                        end
                    end

                    if ( spat | avepat | ipat )

                        if ( option.TwoDPlusTRecon & (REP>1 | PHS>1 | SLC>1 | PAR>1) )

                            if ( REP>1 | PHS>1 ) % use REP or PHS
                                for slc=1:SLC                               
                                    if ( REP > 1 )
                                        for phs=1:PHS

                                            if ( ipat )
                                                kspaceForRecon = kspace_ipat_withRef(:,:,:,acq, slc, 1, eco, phs, :, set);
                                                refForRecon = ref_ipat(:,:,:,acq, slc, 1, eco, phs, :, set);
                                            else                                            
                                                kspaceForRecon = kspace(:,:,:,acq, slc, 1, eco, phs, :, set);
                                                if ( avepat ) 
                                                    refForRecon = ref(:,:,:,1, slc, 1, eco, phs, :, set);
                                                else
                                                    refForRecon = ref(:,:,:,acq, slc, 1, eco, phs, :, set);
                                                end
                                            end

                                            kspaceForRecon = squeeze(kspaceForRecon);
                                            refForRecon = squeeze(refForRecon);

                                            [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS(kspaceForRecon, refForRecon, option);

                                            dstCHA = size(kspaceLinear, 3);

                                            meanKSpace = mean(kspaceLinear, 4);
                                            senMap = CoilSensitivityMapEstimation(ifft2c(meanKSpace), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
                                            imLinear = SensitivityCoilCombination(ifft2c(kspaceLinear), senMap);
                                            imVICPAAS = SensitivityCoilCombination(ifft2c(kspaceVICPAAS), senMap);

                                            for rep=1:REP
                                                fullKSpaceLinear(:,:,1:dstCHA,acq, slc, 1, eco, phs, rep, set) = kspaceLinear(:,:,:,rep);
                                                unwarppedImLinear(:,:,acq, slc, 1, eco, phs, rep, set) = imLinear(:,:,rep);
                                                fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, 1, eco, phs, rep, set) = kspaceVICPAAS(:,:,:,rep);
                                                unwarppedImVICPAAS(:,:,acq, slc, 1, eco, phs, rep, set) = imVICPAAS(:,:,rep);
                                            end
                                        end
                                    else
                                        for rep=1:REP
                                            if ( ipat )
                                                kspaceForRecon = kspace_ipat_withRef(:,:,:,acq, slc, 1, eco, :, rep, set);
                                                refForRecon = ref_ipat(:,:,:,acq, slc, 1, eco, :, rep, set);
                                            else
                                                kspaceForRecon = kspace(:,:,:,acq, slc, 1, eco, :, rep, set);
                                                if ( avepat ) 
                                                    refForRecon = ref(:,:,:,1, slc, 1, eco, :, rep, set);
                                                else
                                                    refForRecon = ref(:,:,:,acq, slc, 1, eco, :, rep, set);
                                                end
                                            end

                                            kspaceForRecon = squeeze(kspaceForRecon);
                                            refForRecon = squeeze(refForRecon);

                                            [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS(kspaceForRecon, refForRecon, option);

                                            dstCHA = size(kspaceLinear, 3);

                                            meanKSpace = mean(kspaceLinear, 4);
                                            senMap = CoilSensitivityMapEstimation(ifft2c(meanKSpace), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
                                            imLinear = SensitivityCoilCombination(ifft2c(kspaceLinear), senMap);
                                            imVICPAAS = SensitivityCoilCombination(ifft2c(kspaceVICPAAS), senMap);

                                            for phs=1:PHS
                                                fullKSpaceLinear(:,:,1:dstCHA,acq, slc, 1, eco, phs, rep, set) = kspaceLinear(:,:,:,phs);
                                                unwarppedImLinear(:,:,acq, slc, 1, eco, phs, rep, set) = imLinear(:,:,phs);
                                                fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, 1, eco, phs, rep, set) = kspaceVICPAAS(:,:,:,phs);
                                                unwarppedImVICPAAS(:,:,acq, slc, 1, eco, phs, rep, set) = imVICPAAS(:,:,phs);
                                            end
                                        end
                                    end
                                end
                            else
                                % use SLC or PAR
                                if ( SLC>1 & PAR==1 )
                                    if ( ipat )
                                        kspaceForRecon = kspace_ipat_withRef(:,:,:,acq, :, 1, eco, 1, 1, set);
                                        refForRecon = ref_ipat(:,:,:,acq, :, 1, eco, 1, 1, set);
                                    else
                                        kspaceForRecon = kspace(:,:,:,acq, :, 1, eco, 1, 1, set);
                                        if ( avepat )
                                            refForRecon = ref(:,:,:,1, :, 1, eco, 1, 1, set);
                                        else
                                            refForRecon = ref(:,:,:,acq, :, 1, eco, 1, 1, set);
                                        end
                                    end

                                    kspaceForRecon = squeeze(kspaceForRecon);
                                    refForRecon = squeeze(refForRecon);

                                    % across SLC, better turn off the temporal reg
                                    option.temporalScalingFactorSPIRiT = 0;   
                                    option.dynamicKernel = 1;
                                    [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS(kspaceForRecon, refForRecon, option);

                                    dstCHA = size(kspaceLinear, 3);

                                    for slc=1:SLC                                                      
                                        fullKSpaceLinear(:,:,1:dstCHA,acq, slc, 1, eco, 1, 1, set) = kspaceLinear(:,:,:,1,slc);
                                        unwarppedImLinear(:,:,acq, slc, 1, eco, 1, 1, set) = unwrappedImVICPAASLinear(:,:,1,slc);

                                        fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, 1, eco, 1, 1, set) = kspaceVICPAAS(:,:,:,1,slc);
                                        unwarppedImVICPAAS(:,:,acq, slc, 1, eco, 1, 1, set) = unwrappedImVICPAAS(:,:,1,slc);
                                    end
                                end

                                if ( SLC==1 & PAR>1 )
                                    if ( ipat )
                                        kspaceForRecon = kspace_ipat_withRef(:,:,:,acq, 1,:, eco, 1, 1, set);
                                        refForRecon = ref_ipat(:,:,:,acq, 1, :, eco, 1, 1, set);
                                    else
                                        kspaceForRecon = kspace(:,:,:,acq, 1, :, eco, 1, 1, set);
                                        if ( avepat )
                                            refForRecon = ref(:,:,:,1, 1, :, eco, 1, 1, set);
                                        else
                                            refForRecon = ref(:,:,:,acq, 1, :, eco, 1, 1, set);
                                        end
                                    end

                                    % COL LIN CHA PAR
                                    kspaceForRecon = squeeze(kspaceForRecon);
                                    refForRecon = squeeze(refForRecon);

%                                     kspaceForRecon = reshape(kspaceForRecon, [size(kspaceForRecon,1) size(kspaceForRecon,2) size(kspaceForRecon,3) 1 size(kspaceForRecon, 6)]);
%                                     refForRecon = reshape(refForRecon, [size(refForRecon,1) size(refForRecon,2) size(refForRecon,3) 1 size(refForRecon, 6)]);
                                    
                                    option.dynamicKernel = 1;
                                    [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS(kspaceForRecon, refForRecon, option);

                                    dstCHA = size(kspaceLinear, 3);

                                    for par=1:PAR                                                      
                                        fullKSpaceLinear(:,:,1:dstCHA,acq, 1, par, eco, 1, 1, set) = kspaceLinear(:,:,:,par);
                                        unwarppedImLinear(:,:,acq, 1, par, eco, 1, 1, set) = unwrappedImVICPAASLinear(:,:,par);
                                        fullKSpaceVICPAAS(:,:,1:dstCHA,acq, 1, par, eco, 1, 1, set) = kspaceVICPAAS(:,:,:,par);
                                        unwarppedImVICPAAS(:,:,acq, 1, par, eco, 1, 1, set) = unwrappedImVICPAAS(:,:,par);
                                    end
                                end

                                if ( SLC>1 & PAR>1 )

                                    % go through every SLC, run recon on PAR

                                    for slc=1:SLC
                                        disp(['slc - ' num2str(slc)]);
                                        if ( ipat )
                                            kspaceForRecon = kspace_ipat_withRef(:,:,:,acq, slc,:, eco, 1, 1, set);
                                            refForRecon = ref_ipat(:,:,:,acq, slc, :, eco, 1, 1, set);
                                        else
                                            kspaceForRecon = kspace(:,:,:,acq, slc, :, eco, 1, 1, set);
                                            if ( avepat )
                                                refForRecon = ref(:,:,:,1, slc, :, eco, 1, 1, set);
                                            else
                                                refForRecon = ref(:,:,:,acq, slc, :, eco, 1, 1, set);
                                            end
                                        end

                                        kspaceForRecon = squeeze(kspaceForRecon);
                                        refForRecon = squeeze(refForRecon);
                                        
%                                         kspaceForRecon = reshape(kspaceForRecon, [size(kspaceForRecon,1) size(kspaceForRecon,2) size(kspaceForRecon,3) 1 size(kspaceForRecon, 6)]);
%                                         refForRecon = reshape(refForRecon, [size(refForRecon,1) size(refForRecon,2) size(refForRecon,3) 1 size(refForRecon, 6)]);
                                        
                                        option.dynamicKernel = 1;
                                        [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS(kspaceForRecon, refForRecon, option);

                                        dstCHA = size(kspaceLinear, 3);

                                        % since the SLC is used, a coil map is need for a slice
                                        senMap = CoilSensitivityMapEstimation(ifft2c(kspaceLinear), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
                                        for par=1:PAR                                                      
                                            %senMap = CoilSensitivity_Souheil_Walsh(ifft2c(kspaceLinear(:,:,:,par)), option.kSizeEigenVectorCoilSensitivity, 0);
                                            imLinear = SensitivityCoilCombination(ifft2c(kspaceLinear(:,:,:,par)), senMap(:,:,:,par));
                                            imVICPAAS = SensitivityCoilCombination(ifft2c(kspaceVICPAAS(:,:,:,par)), senMap(:,:,:,par));

                                            fullKSpaceLinear(:,:,1:dstCHA,acq, slc, par, eco, 1, 1, set) = kspaceLinear(:,:,:,par);
                                            unwarppedImLinear(:,:,acq, slc, par, eco, 1, 1, set) = imLinear;
                                            fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, par, eco, 1, 1, set) = kspaceVICPAAS(:,:,:,par);
                                            unwarppedImVICPAAS(:,:,acq, slc, par, eco, 1, 1, set) = imVICPAAS;
                                        end
                                    end
                                end

                                % plotComplexImageArray(squeeze(unwarppedImVICPAAS(:,:,acq, slc, :, eco, 1, 1, set)), [1 1 1], LIN, 1204, 2048, 0.1, 1);
                                % plotComplexImageArray(squeeze(unwarppedImLinear(:,:,acq, slc, :, eco, 1, 1, set)), [1 1 1], LIN, 1204, 2048, 0.1, 1);
                            end
                        else
                            for slc=1:SLC
                                for rep=1:REP
                                    for phs=1:PHS
                                        if ( ipat )
                                            kspaceForRecon = kspace_ipat_withRef(:,:,:,acq, slc, 1, eco, phs, rep, set);
                                            refForRecon = ref_ipat(:,:,:,acq, slc, 1, eco, phs, rep, set);
                                        else
                                            kspaceForRecon = kspace(:,:,:,acq, slc, 1, eco, phs, rep, set);
                                            if ( avepat ) 
                                                refForRecon = ref(:,:,:,1, slc, 1, eco, phs, rep, set);
                                            else
                                                refForRecon = ref(:,:,:,acq, slc, 1, eco, phs, rep, set);
                                            end
                                        end

                                        [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS(kspaceForRecon, refForRecon, option);

                                        dstCHA = size(kspaceLinear, 3);

                                        meanKSpace = mean(kspaceLinear, 4);
                                        senMap = CoilSensitivityMapEstimation(ifft2c(meanKSpace), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
                                        imLinear = SensitivityCoilCombination(ifft2c(kspaceLinear), senMap);
                                        imVICPAAS = SensitivityCoilCombination(ifft2c(kspaceVICPAAS), senMap);

                                        fullKSpaceLinear(:,:,1:dstCHA,acq, slc, 1, eco, phs, rep, set) = kspaceLinear;
                                        unwarppedImLinear(:,:,acq, slc, 1, eco, phs, rep, set) = imLinear;
                                        fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, 1, eco, phs, rep, set) = kspaceVICPAAS;
                                        unwarppedImVICPAAS(:,:,acq, slc, 1, eco, phs, rep, set) = imVICPAAS;
                                    end
                                end
                            end
                        end

                    end
                end
            end
        end

        unwarppedImVICPAAS = performKSpaceFilter2D(fft2c(unwarppedImVICPAAS), filterCOL, filterLIN);
        unwarppedImVICPAAS = ifft2c(unwarppedImVICPAAS);
        
        % perfrom zero-padding resize
        unwarppedImLinearResized = Zero_Padding_Resize_NoFiltering(unwarppedImLinear, COL, newNpe);
        unwarppedImVICPAASResized = Zero_Padding_Resize_NoFiltering(unwarppedImVICPAAS, COL, newNpe);
    end
else
       
    % 3D recon
    if ( performGRAPPA & ~isRandomSampling )
        
        LIN = size(kspace, 2);
        
        fullKSpaceGrappa = kspace;              
        unwarppedImGrappa = zeros(COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET);
        gFactorGrappa = zeros(COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET);       
        
        for set=1:SET
            for eco=1:ECO
                for slc=1:SLC
                    for acq=1:ACQ

                        disp(['acq - slc - eco - set : ' num2str(acq) ' ' num2str(slc) ' ' num2str(eco) ' ' num2str(set)]);

                        if ( tpat )

                            if ( tpat_useREP )
                                kspaceForRecon = kspace(:,:,:,acq, slc, :, eco, 1, :, set);
                            end

                            if ( tpat_usePHS )
                                kspaceForRecon = kspace(:,:,:,acq, slc, :, eco, :, 1, set);
                            end

                            kspaceForRecon = squeeze(kspaceForRecon);
                            [im, fullkspace, sensitivityMap, gFactor] = TGRAPPA3D_AverageAll_SrcDstChannels_SNRUnit(kspaceForRecon, FOV_reduction_factor, FOV_reduction_factor3D, option);

                            meanKSpace = mean(fullkspace, 5); % COL LIN CHA PAR N
                            meanKSpace = permute(meanKSpace, [1 2 4 3]);% COL LIN PAR CHA N
                            meanKSpace = ifft3c(meanKSpace);
                            meanKSpace = permute(meanKSpace, [1 2 4 3]);% COL LIN CHA PAR N
                            senMap = CoilSensitivityMapEstimation(ifft3c(meanKSpace), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
                            
                            complexIm = permute(fullkspace, [1 2 4 3 5]);% COL LIN PAR CHA N
                            complexIm = ifft3c(complexIm);
                            complexIm = permute(complexIm, [1 2 4 3 5]);% COL LIN CHA PAR N
                            im = SensitivityCoilCombination(complexIm, senMap);

                            dstCha = size(fullkspace, 3);

                            if ( tpat_useREP )
                                for rep=1:REP
                                    fullKSpaceGrappa(:,:,1:dstCha,acq, slc, :, eco, 1, rep, set) = fullkspace(:,:,:,:,rep);
                                    unwarppedImGrappa(:,:,acq, slc, :, eco, 1, rep, set) = im(:,:,:,rep);
                                    gFactorGrappa(:,:,acq, slc, :, eco, 1, rep, set) = gFactor;
                                end
                            end

                            if ( tpat_usePHS )
                                for phs=1:PHS
                                    fullKSpaceGrappa(:,:,1:dstCha,acq, slc, :, eco, phs, 1, set) = fullkspace(:,:,:,:,phs);
                                    unwarppedImGrappa(:,:,acq, slc, :, eco, phs, 1, set) = im(:,:,:,phs);
                                    gFactorGrappa(:,:,acq, slc, :, eco, phs, 1, set) = gFactor;
                                end
                            end
                        end

                        if ( spat | avepat )                        
                            for rep=1:REP
                                for phs=1:PHS
                                    kspaceForRecon = kspace(:,:,:,acq, slc, :, eco, phs, rep, set);

                                    setRef = set;
                                    if ( refSET < set )
                                        setRef = refSET;
                                    end

                                    if (avepat) 
                                        refForRecon = ref(:,:,:,1, slc, :, eco, phs, rep, setRef);
                                    else
                                        if ( refREP < rep )
                                            refForRecon = ref(:,:,:,acq, slc, :, eco, phs, refREP, setRef);
                                        else
                                            refForRecon = ref(:,:,:,acq, slc, :, eco, phs, rep, setRef);
                                        end
                                    end

                                    if ( refREP > 1 )                                        
                                        [im, fullkspace, sensitivityMap, gFactor] = TGRAPPA3D_SeperateRef_SrcDstChannels_SNRUnit(kspaceForRecon, refForRecon, FOV_reduction_factor, FOV_reduction_factor3D, option);
                                        dstCha = size(fullkspace, 3);
                                        fullKSpaceGrappa(:,:,1:dstCha,acq, slc, :, eco, phs, rep, set) = fullkspace;

%                                             senMap = CoilSensitivityMapEstimation(ifft2c(fullkspace), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
%                                             im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);
                                        unwarppedImGrappa(:,:,acq, slc, :, eco, phs, rep, set) = im;
                                        gFactorGrappa(:,:,acq, slc, :, eco, phs, rep, set) = gFactor;
                                    else
                                        kspaceForRecon = kspace(:,:,:,acq, slc, :, eco, phs, :, set);
                                        kspaceForRecon = squeeze(kspaceForRecon);

                                        [im, fullkspace, sensitivityMap, gFactor] = TGRAPPA3D_SeperateRef_SrcDstChannels_SNRUnit(kspaceForRecon, refForRecon, FOV_reduction_factor, FOV_reduction_factor3D, option);

                                        dstCha = size(fullkspace, 3);
                                        fullKSpaceGrappa(:,:,1:dstCha,acq, slc, :, eco, phs, 1:size(fullkspace,4), set) = fullkspace;

%                                             senMap = CoilSensitivityMapEstimation(ifft2c(mean(fullkspace, 4)), option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
%                                             im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);
                                        unwarppedImGrappa(:,:,acq, slc, :, eco, phs, 1:size(im,3), set) = im;
                                        gFactorGrappa(:,:,acq, slc, :, eco, phs, rep, set) = gFactor;
                                    end
                                end

                                if ( refREP == 1 )
                                    break;
                                end
                            end
                        end

                        if ( ipat )
                            disp('ipat ... ');
                            for rep=1:REP
                                for phs=1:PHS
                                    kspaceForRecon = squeeze(kspace(:,:,:,acq, slc, :, eco, phs, rep, set));
                                    refForRecon = squeeze(ref(:,:,:,acq, slc, :, eco, phs, rep, set));

                                    [im, fullkspace, sensitivityMap, gFactor, E_0, V] = TGRAPPA3D_InplaceRef_SrcDstChannels_SNRUnit(kspaceForRecon, refForRecon, FOV_reduction_factor, FOV_reduction_factor3D, option);

                                    dstCha = size(fullkspace, 3);
                                    fullKSpaceGrappa(:,:,1:dstCha,acq, slc, :, eco, phs, rep, set) = fullkspace;

                                    complexIm = permute(fullkspace, [1 2 4 3]);
                                    complexIm = ifft3c(complexIm);
                                    complexIm = permute(complexIm, [1 2 4 3]);
                                    
                                    senMap = CoilSensitivityMapEstimation(complexIm, option.csmMethod, option.spatialSmoothingKernel, option.kSizeEigenVectorCoilSensitivity, option.percentageEigenVectorCoilSensitivity, option.kSizeEigenVectorCoilSensitivity);
                                    
                                    im = SensitivityCoilCombination(complexIm, senMap);                                                                               

                                    unwarppedImGrappa(:,:,acq, slc, :, eco, phs, rep, set) = im;
                                    gFactorGrappa(:,:,acq, slc, :, eco, phs, rep, set) = gFactor;
                                end
                            end
                        end
                    end
                end
            end
        end

        unwarppedImGrappa = performKSpaceFilter3D(fft3c_ImageDimensions(unwarppedImGrappa), filterCOL, filterLIN, filterPAR);
        unwarppedImGrappa = ifft3c_ImageDimensions(unwarppedImGrappa);

        unwarppedImGrappaResized = Zero_Padding_Resize_NoFiltering3D(unwarppedImGrappa, COL, newNpe, newNpar);

        clear unwarppedImGrappa
    end
    
    if ( performVICPAAS )
        
        option.temporalScalingFactorSPIRiT = 2;
        % option.wavWeightSPIRiT = option.wavWeightSPIRiT*2;
        option.dstChaThres  = 0.002;
        option.maxIterSPIRiT = 70;
        option.stopThresSPIRiT = 1e-5;
    
        fullKSpaceLinear = kspace;
        unwarppedImLinear = zeros(COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET);

        fullKSpaceVICPAAS = kspace;
        unwarppedImVICPAAS = unwarppedImLinear;

        for set=1:SET
            for eco=1:ECO
                for acq=1:ACQ

                    if ( spat | avepat | ipat )                    
                        for phs=1:PHS
                            for rep=1:REP
                                for slc=1:SLC
                                    disp(['slc - ' num2str(slc)]);
                                    if ( ipat )
                                        kspaceForRecon = kspace_ipat_withRef(:,:,:,acq, slc,:, eco, phs, rep, set);
                                        refForRecon = ref_ipat(:,:,:,acq, slc, :, eco, phs, rep, set);
                                    else
                                        kspaceForRecon = kspace(:,:,:,acq, slc, :, eco, phs, rep, set);

                                        if ( avepat )
                                            refForRecon = ref(:,:,:,1, slc, :, eco, phs, rep, set);
                                        else
                                            if ( set > size(ref, 10) )
                                                refForRecon = ref(:,:,:,acq, slc, :, eco, phs, rep, end);
                                            else
                                                refForRecon = ref(:,:,:,acq, slc, :, eco, phs, rep, set);
                                            end
                                        end
                                    end

                                    kspaceForRecon = squeeze(kspaceForRecon); % plotKSpaceSamplingPattern(squeeze(kspaceForRecon(end/2,:,1,:)));
                                    refForRecon = squeeze(refForRecon); % plotKSpaceSamplingPattern(squeeze(refForRecon(end/2,:,1,:)));

                                    [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS3D(kspaceForRecon, refForRecon, option, rangeUsedFE);

                                    dstCHA = size(kspaceLinear, 3);

                                    fullKSpaceLinear(:,:,1:dstCHA,acq, slc, :, eco, phs, rep, set) = kspaceLinear;
                                    unwarppedImLinear(:,:,acq, slc, :, eco, phs, rep, set) = unwrappedImVICPAASLinear;
                                    fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, :, eco, phs, rep, set) = kspaceVICPAAS;
                                    unwarppedImVICPAAS(:,:,acq, slc, :, eco, phs, rep, set) = unwrappedImVICPAAS;
                                end                    
                            end
                        end
                    end

                    if ( tpat )                    
                        for slc=1:SLC                        
                            if ( tpat_useREP )
                                kspaceForRecon = kspace(:,:,:,acq, slc, :, eco, 1, :, set);
                                ref = mean(kspaceForRecon, 9);
                            end

                            if ( tpat_usePHS )
                                kspaceForRecon = kspace(:,:,:,acq, slc, :, eco, :, 1, set);
                                ref = mean(kspaceForRecon, 8);
                            end

                            kspaceForRecon = squeeze(kspaceForRecon);    
                            refForRecon = squeeze(ref);

                            if ( tpat_useREP )
                                for rep=1:REP
                                    [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS3D(kspaceForRecon(:,:,:,:,rep), refForRecon, option, rangeUsedFE);
                                    dstCHA = size(kspaceLinear, 3);
                                    fullKSpaceLinear(:,:,1:dstCHA,acq, slc, :, eco, 1, rep, set) = kspaceLinear;
                                    unwarppedImLinear(:,:,acq, slc, 1, eco, :, rep, set) = unwrappedImVICPAASLinear;
                                    fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, :, eco, 1, rep, set) = kspaceVICPAAS;
                                    unwarppedImVICPAAS(:,:,acq, slc, 1, eco, :, rep, set) = unwrappedImVICPAAS;
                                end
                            end

                            if ( tpat_usePHS )
                                for phs=1:PHS
                                    [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS3D(kspaceForRecon(:,:,:,:,phs), refForRecon, option, rangeUsedFE);
                                    dstCHA = size(kspaceLinear, 3);
                                    fullKSpaceLinear(:,:,1:dstCHA,acq, slc, 1, eco, phs, 1, set) = kspaceLinear;
                                    unwarppedImLinear(:,:,acq, slc, 1, eco, phs, 1, set) = unwrappedImVICPAASLinear;
                                    fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, 1, eco, phs, 1, set) = kspaceVICPAAS;
                                    unwarppedImVICPAAS(:,:,acq, slc, 1, eco, phs, 1, set) = unwrappedImVICPAAS;
                                end
                            end
                        end                    
                    end                                                                                
                end
            end
        end    

        unwarppedImVICPAAS = performKSpaceFilter3D(fft3c_ImageDimensions(unwarppedImVICPAAS), filterCOL, filterLIN, filterPAR);
        unwarppedImVICPAAS = ifft3c_ImageDimensions(unwarppedImVICPAAS);

        unwarppedImLinearResized = Zero_Padding_Resize_NoFiltering3D(unwarppedImLinear, COL, newNpe, newNpar);
        unwarppedImVICPAASResized = Zero_Padding_Resize_NoFiltering3D(unwarppedImVICPAAS, COL, newNpe, newNpar);
    end
end

end

unwarppedImGrappaResized = single(unwarppedImGrappaResized);
unwarppedImLinearResized = single(unwarppedImLinearResized);
unwarppedImVICPAASResized = single(unwarppedImVICPAASResized);
fullKSpaceGrappa = single(fullKSpaceGrappa);
fullKSpaceLinear = single(fullKSpaceLinear);
fullKSpaceVICPAAS = single(fullKSpaceVICPAAS);
gFactor = gFactorGrappa;

end

function aZ = addPrePostZeros(centreNo, sampledNo)
    % aZ = 1 : pre zeros
    % aZ = 2 : post zeros
    % aZ = 0 : no zeros
    if ( 2*centreNo == sampledNo )
        aZ = 0;
        return;
    end

    if ( 2*centreNo < sampledNo )
        aZ = 1;
        return;
    end

    if ( 2*centreNo > sampledNo )
        aZ = 2;
        return;
    end        
end