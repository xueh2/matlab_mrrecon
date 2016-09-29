
%% incoming kspace
s = size(kspace)
size(Noise)
sRef = size(ref)

COL = s(1);
LIN = s(2);
CHA = s(3);

%% find the pat mode
RefScanMode = protocol.sPat.ucRefScanMode;
ind = find(RefScanMode=='0');
RefScanMode = RefScanMode(ind:end)

FOV_reduction_factor = protocol.sPat.lAccelFactPE
FOV_reduction_factor3D = protocol.sPat.lAccelFact3D
numRefLines = protocol.sPat.lRefLinesPE

ipat = 0;
spat = 0;
tpat = 0;
avepat = 0;

if ( ~isempty(strfind(RefScanMode, '0x4')) & isempty(strfind(RefScanMode, '0x40')) )
    ipat = 1;
end

if ( ~isempty(strfind(RefScanMode, '0x2')) )
    spat = 1;
end

if ( ~isempty(strfind(RefScanMode, '0x20')) | ~isempty(strfind(RefScanMode, '0x40')) )
    tpat = 1;
end

if ( ~isempty(strfind(RefScanMode, '0x10')) )
    avepat = 1;
end

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

if ( SEG > 1 )
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

    if ( refSEG > 1 )
        ref2 = sum(ref, 11);
        ref = ref2;
        sRef = size(ref)
    end

    maxDimRef = numel(sRef);
end

%% remove read-out oversampling
kspace = performDownSampleFE(kspace);
COL = size(kspace, 1);
size(kspace)

ref = performDownSampleFE(ref);
refCOL = size(ref, 1);
size(ref)

%% noise prewhitening
disp('performing noise prewhitening...') 

Noise = squeeze(Noise);
Noise = permute(Noise, [2 1]);
noisePrewhiteningMatrix = calculateNoisePrewhitener(Noise, rx_dwelltime_data);

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

newNpe = round(peFOV/(feFOV/COL))
if ( LIN > newNpe )
    newNpe = LIN;
end

voxelsize = [feFOV/COL feFOV/COL sliceThickness]

% plotKSpaceArray(kspace);

%% prepare the kspace and ref
%% kspace will be [COL LIN CHA PAR ]

if ( tpat )
    
    tpat_useREP = 0;
    tpat_usePHS = 0;
    
    % REP used as tpat dim
    if ( strfind(RefScanMode, '0x20') )
        tpat_useREP = 1;
    end
    
    % PHS used as tpat dim
    if ( strfind(RefScanMode, '0x40') )
        tpat_usePHS = 1;
    end
    
    option_ori = option;
    option.KernelSizeSPIRiT = option.KernelSizeSPIRiT_TPAT;
    option.OutKernelSPIRiT = option.OutKernelSPIRiT_TPAT;
    option.thresRegSPIRiT = option.thresRegSPIRiT_TPAT;
end

if ( spat )    
    loc = detectSampledLines(ref(:,:,1,1,1,1,1,1,1,1))
    
    stepSize = loc(2:end) - loc(1:end-1);
    
    if ( ~tpat )
        if ( ~isempty(find(stepSize>1)) )
            spat = 0;
            ipat = 1;
        else    
            ref= ref(:,loc(1):loc(end), :,:,:,:,:,:,:,:);
        end
    end
end

if ( ipat )          
    kspace_ipat_withRef = kspace;
    kspace_ipat_withRef(:,1:refLIN,:,:,:,:,:,:,:,:) = ref + kspace(:,1:refLIN,:,:,:,:,:,:,:,:);

    ref_ipat = ref + kspace(:,1:refLIN,:,:,:,:,:,:,:,:);
    loc = detectSampledLines(ref_ipat(:,:,1,1,1,1,1,1,1,1))
    diff = loc([2:end 1]) - loc;
    ind = find(diff==1);
    startLine = loc(ind(1));
    endLine = loc(ind(end))+1;
    ref_ipat= ref_ipat(:,startLine:endLine, :,:,:,:,:,:,:,:);
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

%% 2D or 2D + T recon
if ( FOV_reduction_factor3D == 1 ) 
    
    if ( PAR > 1 )
        kspace = ifftc(kspace,6);
        ref = ifftc(ref,6);
        
        if ( ipat )
            kspace_ipat_withRef = ifftc(kspace_ipat_withRef,6);
            ref_ipat = ifftc(ref_ipat,6);
        end
    end
    
    %% linear
    if ( performGRAPPA )
        fullKSpaceGrappa = kspace;
        unwarppedImGrappa = zeros(COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET);

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
                                senMap = CoilSensitivity_Souheil_Walsh(ifft2c(meanKSpace), option.kSizeEigenVectorCoilSensitivity, 0);
                                im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);

                                dstCha = size(fullkspace, 3);
                                
                                if ( tpat_useREP )
                                    for rep=1:REP
                                        fullKSpaceGrappa(:,:,1:dstCha,acq, slc, par, eco, 1, rep, set) = fullkspace(:,:,:,rep);
                                        unwarppedImGrappa(:,:,acq, slc, par, eco, 1, rep, set) = im(:,:,rep);
                                    end
                                end

                                if ( tpat_usePHS )
                                    for phs=1:PHS
                                        fullKSpaceGrappa(:,:,1:dstCha,acq, slc, par, eco, phs, 1, set) = fullkspace(:,:,:,phs);
                                        unwarppedImGrappa(:,:,acq, slc, par, eco, phs, 1, set) = im(:,:,phs);
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
                                            refForRecon = ref(:,:,:,acq, slc, par, eco, phs, rep, set);
                                        end

                                        [im, fullkspace, sensitivityMap, gFactor] = TGRAPPA_SeperateRef_SrcDstChannels_SNRUnit(kspaceForRecon, refForRecon, FOV_reduction_factor, option);
                                        dstCha = size(fullkspace, 3);
                                        fullKSpaceGrappa(:,:,1:dstCha,acq, slc, par, eco, phs, rep, set) = fullkspace;

                                        senMap = CoilSensitivity_Souheil_Walsh(ifft2c(fullkspace), option.kSizeEigenVectorCoilSensitivity, 0);
                                        im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);
                                        unwarppedImGrappa(:,:,acq, slc, par, eco, phs, rep, set) = im;
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
                                        dstCha = size(fullkspace, 3);
                                        fullKSpaceGrappa(:,:,1:dstCha,acq, slc, par, eco, phs, rep, set) = fullkspace;

                                        senMap = CoilSensitivity_Souheil_Walsh(ifft2c(fullkspace), option.kSizeEigenVectorCoilSensitivity, 0);
                                        im = SensitivityCoilCombination(ifft2c(fullkspace), senMap);
                                        unwarppedImGrappa(:,:,acq, slc, par, eco, phs, rep, set) = im;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        unwarppedImGrappaResized = Zero_Padding_Resize_NoFiltering(unwarppedImGrappa, COL, newNpe);

        clear unwarppedImGrappa
    end
    
    %% non-linear
    if ( performVICPAAS )
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
                                ref = mean(kspaceForRecon, 9);
                            end

                            if ( tpat_usePHS )
                                kspaceForRecon = kspace(:,:,:,acq, slc, 1, eco, :, 1, set);
                                ref = mean(kspaceForRecon, 8);
                            end

                            kspaceForRecon = squeeze(kspaceForRecon);                        

                            if ( option.tspirit & performGRAPPA )
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

                            meanKSpace = mean(kspaceLinear, 4);
                            senMap = CoilSensitivity_Souheil_Walsh(ifft2c(meanKSpace), option.kSizeEigenVectorCoilSensitivity, 0);
                            imLinear = SensitivityCoilCombination(ifft2c(kspaceLinear), senMap);
                            imVICPAAS = SensitivityCoilCombination(ifft2c(kspaceVICPAAS), senMap);

                            if ( tpat_useREP )
                                for rep=1:REP
                                    fullKSpaceLinear(:,:,1:dstCHA,acq, slc, 1, eco, 1, rep, set) = kspaceLinear(:,:,:,rep);
                                    unwarppedImLinear(:,:,acq, slc, 1, eco, 1, rep, set) = imLinear(:,:,rep);
                                    fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, 1, eco, 1, rep, set) = kspaceVICPAAS(:,:,:,rep);
                                    unwarppedImVICPAAS(:,:,acq, slc, 1, eco, 1, rep, set) = imVICPAAS(:,:,rep);
                                end
                            end

                            if ( tpat_usePHS )
                                for phs=1:PHS
                                    fullKSpaceLinear(:,:,1:dstCHA,acq, slc, 1, eco, phs, 1, set) = kspaceLinear(:,:,:,phs);
                                    unwarppedImLinear(:,:,acq, slc, 1, eco, phs, 1, set) = imLinear(:,:,phs);
                                    fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, 1, eco, phs, 1, set) = kspaceVICPAAS(:,:,:,phs);
                                    unwarppedImVICPAAS(:,:,acq, slc, 1, eco, phs, 1, set) = imVICPAAS(:,:,phs);
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
                                            senMap = CoilSensitivity_Souheil_Walsh(ifft2c(meanKSpace), option.kSizeEigenVectorCoilSensitivity, 0);
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
                                            senMap = CoilSensitivity_Souheil_Walsh(ifft2c(meanKSpace), option.kSizeEigenVectorCoilSensitivity, 0);
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
                                    [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS(kspaceForRecon, refForRecon, option);

                                    dstCHA = size(kspaceLinear, 3);

                                    % since the SLC is used, a coil map is need for a slice
                                    senMap = CoilSensitivity_Souheil_Walsh(ifft2c(kspaceLinear), option.kSizeEigenVectorCoilSensitivity, 0);
                                    for slc=1:SLC                                                      
                                        %senMap = CoilSensitivity_Souheil_Walsh(ifft2c(kspaceLinear(:,:,:,slc)), option.kSizeEigenVectorCoilSensitivity, 0);
                                        imLinear = SensitivityCoilCombination(ifft2c(kspaceLinear(:,:,:,slc)), senMap(:,:,:,slc));
                                        imVICPAAS = SensitivityCoilCombination(ifft2c(kspaceVICPAAS(:,:,:,slc)), senMap(:,:,:,slc));

                                        % fullKSpaceLinear = setKSpaceResult(kspaceLinear(:,:,:,slc), fullKSpaceLinear, acq, slc, 1, eco, 1, 1, set);
                                        fullKSpaceLinear(:,:,1:dstCHA,acq, slc, 1, eco, 1, 1, set) = kspaceLinear(:,:,:,slc);
                                        unwarppedImLinear(:,:,acq, slc, 1, eco, 1, 1, set) = imLinear;

                                        % fullKSpaceVICPAAS = setKSpaceResult(kspaceVICPAAS(:,:,:,slc), fullKSpaceVICPAAS, acq, slc, 1, eco, 1, 1, set);
                                        fullKSpaceVICPAAS(:,:,1:dstCHA,acq, slc, 1, eco, 1, 1, set) = kspaceVICPAAS(:,:,:,slc);
                                        unwarppedImVICPAAS(:,:,acq, slc, 1, eco, 1, 1, set) = imVICPAAS;
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

                                    kspaceForRecon = squeeze(kspaceForRecon);
                                    refForRecon = squeeze(refForRecon);

                                    [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS(kspaceForRecon, refForRecon, option);

                                    dstCHA = size(kspaceLinear, 3);

                                    % since the SLC is used, a coil map is need for a slice
                                    senMap = CoilSensitivity_Souheil_Walsh(ifft2c(kspaceLinear), option.kSizeEigenVectorCoilSensitivity, 0);
                                    for par=1:PAR                                                      
                                        %senMap = CoilSensitivity_Souheil_Walsh(ifft2c(kspaceLinear(:,:,:,par)), option.kSizeEigenVectorCoilSensitivity, 0);
                                        imLinear = SensitivityCoilCombination(ifft2c(kspaceLinear(:,:,:,par)), senMap(:,:,:,par));
                                        imVICPAAS = SensitivityCoilCombination(ifft2c(kspaceVICPAAS(:,:,:,par)), senMap(:,:,:,par));

                                        fullKSpaceLinear(:,:,1:dstCHA,acq, 1, par, eco, 1, 1, set) = kspaceLinear(:,:,:,par);
                                        unwarppedImLinear(:,:,acq, 1, par, eco, 1, 1, set) = imLinear;
                                        fullKSpaceVICPAAS(:,:,1:dstCHA,acq, 1, par, eco, 1, 1, set) = kspaceVICPAAS(:,:,:,par);
                                        unwarppedImVICPAAS(:,:,acq, 1, par, eco, 1, 1, set) = imVICPAAS;
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

                                        [unwrappedImVICPAAS, unwrappedImVICPAASLinear, kspaceLinear, kspaceVICPAAS] = PerformReconVICPAAS(kspaceForRecon, refForRecon, option);

                                        dstCHA = size(kspaceLinear, 3);

                                        % since the SLC is used, a coil map is need for a slice
                                        senMap = CoilSensitivity_Souheil_Walsh(ifft2c(kspaceLinear), option.kSizeEigenVectorCoilSensitivity, 0);
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
                                        senMap = CoilSensitivity_Souheil_Walsh(ifft2c(meanKSpace), option.kSizeEigenVectorCoilSensitivity, 0);
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

        % perfrom zero-padding resize
        unwarppedImLinearResized = Zero_Padding_Resize_NoFiltering(unwarppedImLinear, COL, newNpe);
        unwarppedImVICPAASResized = Zero_Padding_Resize_NoFiltering(unwarppedImVICPAAS, COL, newNpe);
    end
else
    
    % 3D recon
    
end

