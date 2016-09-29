       
function [unwrappedIm, fullkspace, sensitivityMap, unwarppedImCombined, dataSize, voxelsize] = FatWaterRecon_GRAPPA(measDatName, kspace, ref, Noise, csmMethod)
% function [unwrappedIm, fullkspace, sensitivityMap, unwarppedImCombined, dataSize, voxelsize] = FatWaterRecon_GRAPPA(measDatName, kspace, ref, Noise)
%
% this function performs the grappa reconstruction for fat water multi-echo datsets
%
% Inputs:
%   measDatName : meas .dat file name
%   kspace : the undersampled kspace [COL LIN CHA ECO REP SET]
%   ref : the seperate reference kspace [COL LIN CHA SET]; if empty, the TGRAPPA is used
%   Noise : the noise signal for pre-whitening [COL CHA]
%   csmMethod : 'Walsh', 'Jun', or 'Souheil'
%
% Output:
%    unwrappedIm : reconstructed uncombined complex image [COL LIN CHA ECO REP SET]
%    fullkspace : reconstructed full kspace [COL LIN CHA ECO REP SET]
%    sensitivityMap : estimated coil sensitivity [COL LIN CHA ECO REP SET]
%    unwarppedImCombined : reconstructed complex image after coil combination [COL LIN ECO REP SET]
%    dataSize : a vector storing the dimensions [Nfe Npe numOfCoil numOfEcho numOfRep numOfSet]
%    voxelsize : the reconstructed voxel size, the zero-filling has been taken care of
%
%     ***************************************
%     *  Hui Xue  (hui-xue@siemens.com)     *
%     ***************************************

if ( nargin < 5 )
    csmMethod = 'Souheil';
end

[headers,protocol]=read_dat_headers(measDatName);
[feFOV, peFOV, sliceThickness] = findFOVFromConfig(headers.Config)
rx_dwelltime_data = protocol.sRXSPEC.alDwellTime{1};

kspace = performDownSampleFE(kspace);

s = size(ref)
if ( numel(s) == 11 )
    ref = sum(ref, 11);
    ref = squeeze(ref);
end
ref = performDownSampleFE(ref);

S0 = size(kspace)
Nfe = S0(1)
Npe = S0(2)

newNpe = round(peFOV/(feFOV/Nfe))
voxelsize = [feFOV/Nfe feFOV/Nfe sliceThickness];

if ( Npe > newNpe )
    newNpe = Npe;
end

numOfCoil = S0(3)
numOfEcho = S0(4)
numOfRep = 1;
numOfSet = 1;

if ( length(S0) == 4 )
    numOfRep = 1
    numOfSet = S0(4)
end

if ( length(S0) == 5 )
    if (  S0(5) > 2 )
        numOfRep = S0(5)
        numOfSet = 1
    else
        numOfRep = 1
        numOfSet = S0(5)
    end
end

if ( length(S0) == 6 )
    numOfRep = S0(5)
    numOfSet = S0(6)
end

dataSize = [Nfe Npe numOfCoil numOfEcho numOfRep numOfSet];
kspace = reshape(kspace, dataSize);
size(kspace)

%% recon         
sampling_location = detectSampledLinesDynamic(kspace);    
FOV_reduction_factor = sampling_location(2,1)-sampling_location(1,1);
    
sampling_location_ref = detectSampledLinesDynamic(ref)
ref = ref(:, sampling_location_ref, :);

% -----------------------------------------
% noise
disp('performing noise prewhitening...') 

noisePK = permute(Noise, [2 1]);
size(noisePK)

noisePrewhiteningMatrix = calculateNoisePrewhitener(noisePK, rx_dwelltime_data);

kspace2 = permute(kspace, [3 1 2 4 5 6]);
kspace2 = applyNoisePrewhitener(kspace2, noisePrewhiteningMatrix);
kspace2 = permute(kspace2, [2 3 1 4 5 6]);
size(kspace2)
kspace = kspace2;
size(kspace)
clear kspace2

if ( ~isempty(ref) )
    s = size(ref);
    ref2 = permute(ref, [3 1 2 4 5 6]);
    ref2 = applyNoisePrewhitener(ref2, noisePrewhiteningMatrix);
    ref2 = permute(ref2, [2 3 1 4 5 6]);
    size(ref2)
    ref = ref2;
    size(ref)
    clear ref2
end

% ----------------------------
%  grappa part
option.KernelSize = [5 4]
option.thresReg  = 1e-4
option.KernelPattern = [-FOV_reduction_factor:FOV_reduction_factor:2*FOV_reduction_factor];
option.OutPattern = [0:FOV_reduction_factor-1];
option.GrappaOnly = 0;

% common parameters
option.dstCha  = -1;
option.dstChaThres  = -1;

option.rawFilterFE = 'Gaussian';
option.rawFilterFEStrength = 'Weak';
option.rawFilterPE = 'Gaussian';
option.rawFilterPEStrength = 'Weak';

option.acsFilterFE = 'Hanning';
option.acsFilterFEStrength = 'Strong';
option.acsFilterPE = 'Hanning';
option.acsFilterPEStrength = 'Strong';

if ( numOfRep == 1 )
    option.KLTSen = 0;
else
    option.KLTSen = 1;
end

option.numOfModesKeptKLTSen = 3;
option.zeroFilledSize = []
option.spatialSmoothingKernel  = 5;
option.kspaceCenterFE = Nfe/2;

option.csmMethod = csmMethod;

unwrappedIm = zeros([Nfe newNpe numOfCoil numOfEcho numOfRep numOfSet]);
fullkspace = zeros([Nfe newNpe numOfCoil numOfEcho numOfRep numOfSet]);
sensitivityMap = zeros([Nfe newNpe numOfCoil numOfEcho numOfRep numOfSet]);
unwarppedImCombined = zeros([Nfe newNpe numOfEcho numOfRep numOfSet]);
    
for set=1:numOfSet
  
    %% linear recon
    reduced_k_data_Used = zeros([Nfe Npe size(kspace,3) numOfEcho*numOfRep]);

    for pp=1:numOfRep
        reduced_k_data_Used(:,:,:,(pp-1)*numOfEcho+1:pp*numOfEcho) = kspace(:,:,:,:,pp,set);
    end
    
    option2 = option;
    option2.dstCha = -1;
    option2.dstChaThres = -1;
    option2.zeroFilledSize = [Nfe newNpe]
    if ( isempty(ref) )            
        [unwrappedIm_Full, fullkspace_Full, sensitivityMap_Full, gFactor_Full, E_0_Full, V_Full] = TGRAPPA_AverageAll_SrcDstChannels_NotCoilCombine_SNRUnit(reduced_k_data_Used, FOV_reduction_factor, option2);            
    else
        [unwrappedIm_Full, fullkspace_Full, sensitivityMap_Full, gFactor_Full, E_0_Full, V_Full] = TGRAPPA_SeperateRef_SrcDstChannels_NotCoilCombine_SNRUnit(reduced_k_data_Used, ref, FOV_reduction_factor, option2);
    end
       
    for pp=1:numOfRep
        unwrappedIm(:,:,:,:,pp,set) = unwrappedIm_Full(:,:,:,(pp-1)*numOfEcho+1:pp*numOfEcho);
        fullkspace(:,:,:,:,pp,set) = fullkspace_Full(:,:,:,(pp-1)*numOfEcho+1:pp*numOfEcho);
        sensitivityMap(:,:,:,:,pp,set) = sensitivityMap_Full(:,:,:,(pp-1)*numOfEcho+1:pp*numOfEcho);
        unwarppedImCombined(:,:,:,pp,set) = SensitivityCoilCombination(unwrappedIm(:,:,:,:,pp,set), sensitivityMap(:,:,:,:,pp,set));
    end    
end

dataSize = [Nfe newNpe numOfCoil numOfEcho numOfRep numOfSet];

