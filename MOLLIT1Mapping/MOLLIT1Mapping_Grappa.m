
function [unwrappedIm, fullkspace, sensitivityMap, unwrappedImSorted, gFactor, voxelsize, TI, TISorted, AcqTime, AcqTimeSorted] = ...
    MOLLIT1Mapping_Grappa(measDatName, kspace, ref, Noise, csmMethod, ulTimeStamp, ulPMUTimeStamp, sLC, kspaceASCIndex)
% function [unwrappedIm, fullkspace, sensitivityMap, unwarppedImCombined, dataSize, voxelsize] = MOLLIT1Mapping_Grappa(measDatName, kspace, ref, Noise, csmMethod)
%
% this function performs the grappa reconstruction for MOLLI T1 mapping datsets
%
% Inputs:
%   measDatName : meas .dat file name
%   kspace : the undersampled kspace [COL LIN CHA REP]
%   ref : the seperate reference kspace [COL LIN CHA]; if empty, the TGRAPPA is used
%   Noise : the noise signal for pre-whitening [COL CHA]
%   csmMethod : 'Walsh', 'Jun', or 'Souheil'
%
% Output:
%    unwrappedIm : reconstructed complex image [COL LIN REP]
%    fullkspace : reconstructed full kspace [COL LIN CHA REP]
%    sensitivityMap : estimated coil sensitivity [COL LIN CHA REP]
%    unwrappedImSorted : reconstructed complex image after coil combination and sorted by acending order of TIs [COL LIN CHA REP]
%    dataSize : a vector storing the dimensions [Nfe Npe numOfCoil numOfRep]
%    voxelsize : the reconstructed voxel size, the zero-filling has been taken care of
%    TI : inversion time for every REP
%
%     ***************************************
%     *  Hui Xue  (hui-xue@siemens.com)     *
%     ***************************************

if ( nargin < 5 )
    csmMethod = 'Souheil';
end

datfilename = measDatName;
[headers,protocol_header]=read_dat_headers(datfilename);
TI = read_molli_TIs_from_datfile(datfilename)

% find the acqusition time
ind = find(kspaceASCIndex>0);
ind = kspaceASCIndex(ind(:));
ulTimeStampKspace = ulTimeStamp(ind(:));

stepTimeStamp = ulTimeStampKspace(2:end) - ulTimeStampKspace(1:end-1);

ind = find(stepTimeStamp > 50);
if ( numel(ind) ~= numel(TI)-1 )
    error('Number of frames is wrong');
end

AcqTime = zeros(numel(TI), 1);

start = 1;
for f=1:numel(ind)
    ts = ulTimeStampKspace(start:ind(f));
    AcqTime(f) = mean(ts);
    start = ind(f)+1;
end

ts = ulTimeStampKspace(ind(end)+1:end);
AcqTime(numel(TI)) = mean(ts);

AcqTime = AcqTime*2.5; % to ms

% find the pixel size and PE lines
[feFOV, peFOV, sliceThickness] = findFOVFromConfig(headers.Config)
rx_dwelltime_data = protocol_header.sRXSPEC.alDwellTime{1}; % from MeasYaps ascii header

ref = squeeze(ref);

sampledRangeKSpace = detectSampledRangeFE(kspace)
kspace = performDownSampleFE(kspace);

sampledRangeRef = detectSampledRangeFE(ref)
ref = performDownSampleFE(ref);

% if necessary, combine ref with kspace
sampling_location = detectSampledLinesDynamic(kspace);    
sampling_location_ref = detectSampledLines(ref(:,:,:,1));

originalRef = ref;

notFullKSpaceInRef = 0;
if ( sampling_location_ref(2,1)-sampling_location_ref(1,1) > 1 )
    notFullKSpaceInRef = 1;
end

S0 = size(kspace);

Nfe = S0(1);
Npe = S0(2);
numOfCoil = S0(3);
numOfRep = S0(4);

newNpe = round(peFOV/(feFOV/Nfe))
voxelsize = [feFOV/Nfe feFOV/Nfe sliceThickness];

if ( Npe > newNpe )
    newNpe = Npe;
end

%% recon
sampling_location = detectSampledLinesDynamic(kspace);    
FOV_reduction_factor = sampling_location(2,1)-sampling_location(1,1);

% -----------------------------------------
% noise
% -----------------------------------------
disp('performing noise prewhitening...') 

noisePK = permute(Noise, [2 1 3]);
noisePrewhiteningMatrix = calculateNoisePrewhitener(noisePK, rx_dwelltime_data);
reduced_k_data2 = permute(kspace, [3 1 2 4]);
reduced_k_data2 = applyNoisePrewhitener(reduced_k_data2, noisePrewhiteningMatrix);
reduced_k_data2 = permute(reduced_k_data2, [2 3 1 4]);
kspace = reduced_k_data2;
clear reduced_k_data2

ref2 = permute(ref, [3 1 2 4]);
ref2 = applyNoisePrewhitener(ref2, noisePrewhiteningMatrix);
ref2 = permute(ref2, [2 3 1 4]);

ref = ref2;
clear ref2

% ----------------------------
%  grappa part
option.KernelSize = [5 4]
option.thresReg  = 1e-4
option.KernelPattern = [-FOV_reduction_factor:FOV_reduction_factor:2*FOV_reduction_factor];
option.OutPattern = [0:FOV_reduction_factor-1];
option.GrappaOnly = 0;

% v-spirit part
option.KernelSizeSPIRiT = [9 9];
option.OutKernelSPIRiT = [1 1];
option.thresRegSPIRiT = 0.05;
option.maxIterSPIRiT = 70;
option.wavWeightSPIRiT = 0.1;
option.cgSPIRiTlambda = 1e-2;
option.TVWeightSPIRiT = 1e-2;
option.showIterSPIRiT = 0;
option.stopThresSPIRiT = 1e-5;
option.pocsFlagSPIRiT = 0;
option.continuationStep = 10;
option.wavThresRatio = 1;
option.unwarpMethod = 'fft_AllGPU';
option.Itnlim = 5;
option.dataWeightSPIRiT = -1;
option.objTollSPIRiT = 0.1;

option.csmMethod = csmMethod;

% common parameters
option.dstCha  = -1;
option.dstChaThres  = -1;

option.rawFilterFE = 'Gaussian';
option.rawFilterFEStrength = 'Weak';
option.rawFilterPE = 'Gaussian';
option.rawFilterPEStrength = 'Weak';

option.acsFilterFE = 'None';
option.acsFilterFEStrength = 'Strong';
option.acsFilterPE = 'None';
option.acsFilterPEStrength = 'Strong';

option.KLTSen = 1;
option.numOfModesKeptKLTSen = 3;
option.zeroFilledSize = []
option.spatialSmoothingKernel  = 5;
option.kspaceCenterFE = Nfe/2;

option

%% tpat       
kspace2 = zpad(kspace, Nfe, newNpe, size(kspace,3), size(kspace,4));

if ( isempty(ref) )            
    [unwrappedIm, fullkspace, sensitivityMap, gFactor] = TGRAPPA_AverageAll_SrcDstChannels_SNRUnit(kspace2, FOV_reduction_factor, option);
else
    if ( notFullKSpaceInRef )
        [unwrappedIm, fullkspace, sensitivityMap, gFactor, E_0, V] = TGRAPPA_InplaceRef_SrcDstChannels_SNRUnit(kspace2, ref, FOV_reduction_factor, option);
        unwrappedIm = SensitivityCoilCombination(unwrappedIm, sensitivityMap);
    else
        ind = detectSampledLines(ref);
        ref = ref(:,ind, :);
        [unwrappedIm, fullkspace, sensitivityMap, gFactor] = TGRAPPA_SeperateRef_SrcDstChannels_SNRUnit(kspace2, ref, FOV_reduction_factor, option);
    end
end

complexImage = mean(ifft2c(fullkspace), 4);

kSize = 7;
header = CreateFtkHeaderInfo(complexImage, [1 1 1]);
sensitivityMap_final = Matlab_PerformCoilMapEstimation( single(complexImage), header, 'Souheil', 3, kSize, 90);

unwrappedIm = SensitivityCoilCombination(ifft2c(fullkspace), sensitivityMap_final);

% sort the image order
[TISorted, ind] = sort(TI);
unwrappedImSorted = zeros(size(unwrappedIm));
for f=1:numOfRep
    unwrappedImSorted(:,:,f) = unwrappedIm(:,:,ind(f));
end

AcqTimeSorted = AcqTime(ind);
