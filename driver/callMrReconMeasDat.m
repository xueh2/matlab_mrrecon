
% function [unwarppedImGrappaResized, unwarppedImLinearResized, unwarppedImVICPAASResized, fullKSpaceGrappa, fullKSpaceLinear, fullKSpaceVICPAAS, voxelsize] = ...
%             callMrReconMeasDat(datFileName)

% perform the mr recon with vicpaas from a meas .dat file

[folder, name, ext] = fileparts(datFileName);
cd(folder)

dstDir = fullfile(folder, name);
mkdir(dstDir)

h5FileName = fullfile(dstDir, [name '.h5']);

if ( ~isFileExist(h5FileName) )
    command = ['siemens_to_HDF5 ' '"' datFileName '"' ' ' '"' fullfile(dstDir, [name '.h5']) '"']
    dos(command, '-echo');
end

cd(dstDir)
filename = [name '.h5'];

matName = [name '_ori.mat'];

VBorVD = 1;

try
    [headers,protocol]=read_dat_headers(datFileName);
    [feFOV, peFOV, sliceThickness] = findFOVFromConfig(headers.Config)

    sliceThickness = protocol.sSliceArray.asSlice{1}.dThickness;

    save MeasHeader headers protocol

    rx_dwelltime_data = protocol.sRXSPEC.alDwellTime{1}; % from MeasYaps ascii header
    FOV_reduction_factor = protocol.sPat.lAccelFactPE
    FOV_reduction_factor_PAR = protocol.sPat.lAccelFact3D
    protocol.sPat
    
    if ( ~isFileExist(matName) )
        [Data, asc, prot] =Read_RawData(datFileName);
    end
catch
    % if it is VD line
    datfilename = [name '.dat'];
    h5filename = [name '.h5'];

    scanno = 2;
    try        
    [headers, prot] = readSiemensDatHeaderVD(h5filename, scanno);
    catch
        try
        scanno = 1;
        [headers, prot] = readSiemensDatHeaderVD(h5filename, scanno);
        catch
            scanno = 0;
            [headers, prot] = readSiemensDatHeaderVD(h5filename, scanno);
        end
    end

    if ( ~isFileExist(matName) )
        [headers, Data, asc, dataH5] = readSiemensDataVD(h5filename, scanno);
    end
    
    [feFOV, peFOV, sliceThickness] = findFOVFromConfig(headers.Config)

    [MrParcRaidFileHeader, MrParcRaidFileEntry, headers2, protocols]=read_vd_datfile_headers(fullfile(folder, datfilename));   
    sliceThickness = protocols{scanno+1}.sSliceArray.asSlice_{1}.dThickness       
    rx_dwelltime_data = protocols{scanno+1}.sRXSPEC.alDwellTime_{1}
    FOV_reduction_factor = protocols{scanno+1}.sPat.lAccelFactPE
    FOV_reduction_factor_PAR = protocols{scanno+1}.sPat.lAccelFact3D
    protocols{scanno+1}.sPat
    protocol = protocols{scanno+1};
    VBorVD = 2;
    totalAcqTimeInSeconds = 2.5 * double(asc(end).ulTimeStamp - asc(1).ulTimeStamp) * 1e-3;
end

if ( isFileExist(matName) )
    load(matName);
else   
    % [Col Line Cha Acq Slice Partition Echo Phase Rep Set Seg]
%         [Data, asc, prot] =Read_RawData(filename);
%         [kspace, kspaceASCIndex, Noise, ref, refASCIndex, phsCorr, reflect, reflectPhsCorr] = Parse_ICE_Rawdata(Data, asc);
  
    fid = fopen('prot2.txt', 'w');
    fprintf(fid, '%s', char(prot));
    fclose(fid);
     
    [kspace, kspaceASCIndex, Noise, ref, refASCIndex, phsCorr, phsCorrASCIndex, reflect, reflectRef, reflectPhsCorr, other, otherASCIndex, SamplesInScan, KSpaceCentreColumn, MaxKSpaceLineNo, KSpaceCentreLineNo, KSpaceCentrePartitionNo] = Parse_ICE_Rawdata_VD(Data, asc) ;

    if ( size(kspace, 2) > 2*KSpaceCentreLineNo )
        kspace = kspace(:,1:2*KSpaceCentreLineNo,:,:,:,:,:,:,:,:,:);
    end
    
    kspaceOri = kspace;
    refOri = ref;
    phsCorrOri = phsCorr;

    if (0)
        kspace = kspaceOri;
        ref = refOri;
        phsCorr = phsCorrOri;
    end
    
    kspace = PerformRampSamplingOffCenterCorrection(kspace, asc, kspaceASCIndex, reflect, headers, protocol);
    ref = PerformRampSamplingOffCenterCorrection(ref, asc, refASCIndex, reflectRef, headers, protocol);
    phsCorr = PerformRampSamplingOffCenterCorrection(phsCorr, asc, phsCorrASCIndex, reflectPhsCorr, headers, protocol);
    
    kspace = PerformRampSamplingRegridding(kspace, asc, kspaceASCIndex, reflect, headers, protocol);
    ref = PerformRampSamplingRegridding(ref, asc, refASCIndex, reflectRef, headers, protocol);
    phsCorr = PerformRampSamplingRegridding(phsCorr, asc, phsCorrASCIndex, reflectPhsCorr, headers, protocol);
    
    % other = PerformRampSamplingOffCenterCorrection(other, asc, otherASCIndex, headers, protocol);
    
    save([name '_ori.mat'], 'kspace', 'kspaceASCIndex', 'ref', 'refASCIndex', 'Noise', 'phsCorr', 'reflect', 'reflectRef', 'reflectPhsCorr', 'other', 'otherASCIndex', 'prot', 'SamplesInScan', 'KSpaceCentreColumn', 'MaxKSpaceLineNo', 'KSpaceCentreLineNo', 'KSpaceCentrePartitionNo');   

%     kspaceASCIndex = single(kspaceASCIndex);
%     refASCIndex = single(refASCIndex);        
%     save([name '_asc_ori.mat'], 'kspaceASCIndex', 'refASCIndex', 'asc', 'prot');
end

% [Col Line Cha Acq Slice Partition Echo Phase Rep Set Seg]
% 11 dimension array

% kspaceOri = kspace;
% refOri = ref;
% NoiseOri = Noise;
% 
% if ( 0 )
%    kspace = kspaceOri;
%    ref = refOri;
%    Noise = NoiseOri;
% end

%% parameters
if ( ~exist('coilMapMethod') )
    coilMapMethod = 'Souheil';
end

if ( ~exist('dstChaThres') )
    dstChaThres = 0.001;
end

if ( ~exist('performSENSE') )
    performSENSE = 0;
end

if ( ~exist('performGRAPPA') )
    performGRAPPA = 1;
end

if ( ~exist('performVICPAAS') )
    performVICPAAS = 0;
end

if ( VBorVD == 1 )
    rx_dwelltime_noise = bw2rxdwell(130, size(Noise, 1)/2)
else
    rx_dwelltime_noise = 1000*76800.0/size(Noise, 1)/10.0;
end

option = CreateMrReconOption(FOV_reduction_factor, FOV_reduction_factor_PAR, coilMapMethod, dstChaThres, performSENSE, performGRAPPA, performVICPAAS, SamplesInScan, KSpaceCentreColumn, MaxKSpaceLineNo, KSpaceCentreLineNo, KSpaceCentrePartitionNo);

option.PerformPhaseCorrection = 1;
option.PhaseCorrMethod = 'LinearPhaseWithOffset'; % 'PointByPoint', 'LinearPhase', 'LinearPhaseWithOffset'
[unwarppedImSenseResized, unwarppedImGrappaResized, unwarppedImLinear, unwarppedImVICPAA, fullKSpaceGrappa, fullKSpaceLinear, fullKSpaceVICPAAS, voxelsize, gFactorGrappa, gFactorSense] = callMrRecon(kspace, Noise, rx_dwelltime_noise, ref, phsCorr, reflect, reflectRef, reflectPhsCorr, protocol, headers, option, VBorVD);

if ( option.performSENSE )
    % save the results in the order of [Col Line Echo Phase Rep Partition Slice Acq Set]
    if ( ~isempty(unwarppedImSenseResized) )
        unwarppedIm = permute(unwarppedImSenseResized, [1 2 6 7 8 5 4 3 9]);
        unwarppedIm = squeeze(unwarppedIm);
        s = size(unwarppedIm);
        numOfIm = prod(s(3:end));
        unwarppedIm = reshape(unwarppedIm, [s(1) s(2) numOfIm]);
        header = CreateFtkHeaderInfo(unwarppedIm, [voxelsize 1]);
        Matlab_SaveAnalyze(single(abs(unwarppedIm)), header, fullfile(dstDir, 'Sense.hdr'));
        
        save(fullfile(dstDir, 'Sense.mat'), 'unwarppedIm', 'gFactorSense');   
    end
end

if ( option.performGRAPPA )
    % save the results in the order of [Col Line Echo Phase Rep Partition Slice Acq Set]
    if ( ~isempty(unwarppedImGrappaResized) )
        unwarppedIm = permute(unwarppedImGrappaResized, [1 2 6 7 8 5 4 3 9]);
        unwarppedIm = squeeze(unwarppedIm);
        s = size(unwarppedIm);
        numOfIm = prod(s(3:end));
        unwarppedIm = reshape(unwarppedIm, [s(1) s(2) numOfIm]);
        header = CreateFtkHeaderInfo(unwarppedIm, [voxelsize 1]);
        Matlab_SaveAnalyze(single(abs(unwarppedIm)), header, fullfile(dstDir, 'Grappa.hdr'));
        
        save(fullfile(dstDir, 'Grappa.mat'), 'unwarppedIm', 'gFactorGrappa');   
    end
end

if ( option.performVICPAAS )
    unwarppedIm = permute(unwarppedImLinear, [1 2 6 7 8 5 4 3 9]);
    unwarppedIm = squeeze(unwarppedIm);
    s = size(unwarppedIm);
    numOfIm = prod(s(3:end));
    unwarppedIm = reshape(unwarppedIm, [s(1) s(2) numOfIm]);
    header = CreateFtkHeaderInfo(unwarppedIm, [voxelsize 1]);
    Matlab_SaveAnalyze(single(abs(unwarppedIm)), header, fullfile(dstDir, 'VICPAALinear.hdr'));

    save(fullfile(dstDir, 'VICPAALinear.mat'), 'unwarppedIm');   
    
    unwarppedIm = permute(unwarppedImVICPAA, [1 2 6 7 8 5 4 3 9]);
    unwarppedIm = squeeze(unwarppedIm);
    s = size(unwarppedIm);
    numOfIm = prod(s(3:end));
    unwarppedIm = reshape(unwarppedIm, [s(1) s(2) numOfIm]);
    header = CreateFtkHeaderInfo(unwarppedIm, [voxelsize 1]);
    Matlab_SaveAnalyze(single(abs(unwarppedIm)), header, fullfile(dstDir, 'VICPAAS.hdr'));
    
    save(fullfile(dstDir, 'VICPAAS.mat'), 'unwarppedIm');   
end

ref2 = squeeze(ref);
ref2 = ref2(1:128, 1:36, :, 1);
ref2 = single(ref2);

headerSrc = CreateFtkHeaderInfo(ref2, [1 1 1 1]);

fitAcquired = 1;
[kernel, kernelIm, unmixingCoeff, gFactor] = Matlab_PerformVICPAAGKernelCalibrationSrcDst2D(ref2, headerSrc, ref2, headerSrc, 5, [-2 0 2 4], 1, 0, [256 192], 1e-4, 2, fitAcquired);
gFactor1 = gFactor/2;

fitAcquired = 0;
[kernel, kernelIm, unmixingCoeff, gFactor] = Matlab_PerformVICPAAGKernelCalibrationSrcDst2D(ref2, headerSrc, ref2, headerSrc, 5, [-2 0 2 4], 1, 0, [256 192], 1e-4, 2, fitAcquired);
gFactor2 = gFactor/2;

figure; imagescn(cat(3, gFactor1, gFactor2), [0 1.5]);

cd D:\gtuser\gt_windows_setup\ut\VB17A_Gadgetron_testdata\7T\meas_MID709_CV_me_Gadgetron_rep20_bw1002_xRes160_nopat_FID50648

d = ismrmrd.IsmrmrdDataset('meas_MID709_CV_me_Gadgetron_rep20_bw1002_xRes160_nopat_FID50648.h5', '/files/0');
