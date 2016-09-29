
%% ------------------------
% validate 3D recon
% load meas_MID91_t1_fl3d_tra_dynaVIEWS_p2p2_FID2026_ori
% rx_dwelltime_data = 2900;
% 
% load meas_MID5799_T1_VIBE_FS_separate_2x2_FID50187_ori
% rx_dwelltime_data = 3500;

cd D:\data\FromRandall\Kellman\06-3DPSIR-PAT3x2-NotWorking-NeedsLinuxBox\20100413_09h42m49s_7121
load 20100413_09h42m49s_7121_ori
rx_dwelltime_data = 3500

kspace = squeeze(kspace);
ref = squeeze(ref);

kspaceOri = kspace;
refOri = ref;

s = size(kspace)
sRef = size(ref)

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

s = size(kspace);
Nfe = s(1);
Npe = s(2);
CHA = s(3);
Npar = s(4);

sRef = size(ref);
NfeRef = sRef(1);
NpeRef = sRef(2);
CHA = sRef(3);
NparRef = sRef(4);

% coil reduction
[ref2, E_0, V, dstCha] = coilReduction(ref, -1, 1e-3);
kspace2 = applyEigenVector(kspace, V);

ref = ref2;
kspace = kspace2;

kCalib = ref;
kCalib = permute(kCalib, [1 2 4 3]);

plotKSpaceSamplingPattern(kCalib(:,:,1,1));
plotKSpaceSamplingPattern(squeeze(kCalib(1,:,:,1)));

sampledLineLoc = detectSampledLines(kCalib(:,:,1,1))
sampledLineLoc2 = detectSampledLines(squeeze(kCalib(1,:,:,1)))

kCalib = kCalib(:,sampledLineLoc(:),sampledLineLoc2(:),:);
plotKSpaceSamplingPattern(kCalib(:,:,1,1));

size(kCalib)

kSize = [5 5 5];
thresReg = 0.005;
kernel = spiritCalibration3D(kCalib, kSize, thresReg);

GOP = SPIRiT3D(kernel, 'fft', [Nfe, Npe, Npar]);

kspaceUsed = kspace;
kspaceUsed = permute(kspaceUsed, [1 2 4 3]);

plotKSpaceSamplingPattern(kspaceUsed(:,:,1,1));
plotKSpaceSamplingPattern(squeeze(kspaceUsed(1,:,:,1)));

Itnlim = 200;
[res, RESVEC] = cgSPIRiT3D(kspaceUsed, GOP, Itnlim, 1e-6, kspaceUsed);

Im = ifft3c(kspaceUsed);
Img = SoS_Image(Im);
imagescn(real(Img), [], [], [], 3);

Im = ifft3c(kCalib);
Img = SoS_Image(Im);
imagescn(real(Img), [], [], [], 3);

Im = ifft3c(res);
Img = SoS_Image(Im);
imagescn(real(Img), [], [], [], 3);

% try the decoupled recon
kspaceUsed2 = ifftc(kspaceUsed,1);

scale_fctr = norm(kspaceUsed2(:))/(size(kspaceUsed2, 3)*sqrt(CHA));
kspaceUsed2 = kspaceUsed2/scale_fctr;

option.performCenteredFFT = 1;
option.maxIterSPIRiT = 150;
option.stopThresSPIRiT = 1e-6;

option.Itnlim = 10;
option.objTollSPIRiT = 0.01;
option.wavWeightSPIRiT = 0.0025;
option.TVWeightSPIRiT = 0;
option.dataWeightSPIRiT = -1;
option.temporalScalingFactorSPIRiT = 1;
option.TwoDPlusTRecon = 1;
option.UseCoilSensitivity = 1;
option.printInfo = 1;

kspaceUsed3 = permute(kspaceUsed2, [2 3 4 1]);
kernelIm = permute(GOP.KERNEL, [2 3 4 5 1]);
header = CreateFtkHeaderInfo(kspaceUsed3, [1 1 1 1]);
headerKernel = CreateFtkHeaderInfo(kernelIm, [1 1 1 1]);
kspaceLinear = Matlab_PerformSPIRiTLinearRecon2D(single(kspaceUsed3), header, single(kernelIm), headerKernel, single(kspaceUsed3), performCenteredFFT, maxIterSPIRiT, stopThresSPIRiT, printInfo);

kSizeEigenVectorCoilSensitivity = [5 5];
currS_Rep = CoilSensitivity_Souheil_Walsh(ifft2c(kspaceLinear), kSizeEigenVectorCoilSensitivity, 0);
currS_Rep = single(currS_Rep);

res = Matlab_PerformSPIRiTNonLinearRecon2D(single(kspaceUsed3), header, single(kernelIm), headerKernel, ... 
            single(kspaceLinear), single(currS_Rep), option.performCenteredFFT, option.Itnlim, option.objTollSPIRiT, option.wavWeightSPIRiT, ...
            1, option.TVWeightSPIRiT, option.dataWeightSPIRiT, option.temporalScalingFactorSPIRiT, option.TwoDPlusTRecon, 0, ...
            option.UseCoilSensitivity, option.printInfo); 


ImLinear = SensitivityCoilCombination(ifft2c(kspaceLinear), currS_Rep);
plotComplexImageArray(ImLinear);

ImNonLinear = SensitivityCoilCombination(ifft2c(res), currS_Rep);
plotComplexImageArray(ImNonLinear);

plotKSpaceArray(kspaceLinear);
plotKSpaceArray(res);

plotKSpace(kspaceLinear(:,:,:, 25));
plotKSpace(res(:,:,:, 25));

ImLinear = permute(ImLinear, [3 1 2]);
plotComplexImageArray(ImLinear);

ImNonLinear = permute(ImNonLinear, [3 1 2]);

imagescn(abs(ImLinear), [], [], [], 3);
figure; imagescn(abs(ImNonLinear), [0 2], [], [], 3);

for f=1:Nfe
    f
    aKernelIm = GOP.KERNEL(f,:, :, :, :);
    aKernelIm = squeeze(aKernelIm);

    kspaceFrame = kspaceUsed2(f,:,:,:);
    kspaceFrame = squeeze(kspaceFrame);

%     GOP2 = SPIRiT_KernelIm(aKernelIm, 'fft_AllGPU');
%     [kspaceLinear, RESVEC] = cgSPIRiT_GPU(kspaceFrame, GOP2, 200, 1e-6, kspaceFrame);
% 
%     res(f,:,:,:) = kspaceLinear;   
    
    header = CreateFtkHeaderInfo(kspaceFrame, [1 1 1 1]);
    headerKernel = CreateFtkHeaderInfo(aKernelIm, [1 1 1 1]);
    kspaceLinear2 = Matlab_PerformSPIRiTLinearRecon2D(single(kspaceFrame), header, single(aKernelIm), headerKernel, single(kspaceFrame), performCenteredFFT, maxIterSPIRiT, stopThresSPIRiT, printInfo);

    res(f,:,:,:) = kspaceLinear2;
end

res = ifftc(res,2);
res = ifftc(res,3);

Img = SoS_Image(res);
imagescn(real(Img), [], [], [], 3);


plotKSpace(kspaceFrame);

plotKSpace(kspaceLinear);

header = CreateFtkHeaderInfo(kspaceDst, [1 1 1 1]);
headerKernel = CreateFtkHeaderInfo(G_I_Im_REP, [1 1 1 1]);
kspaceLinear = Matlab_PerformSPIRiTLinearRecon2D(single(kspaceDst), header, single(G_I_Im_REP), headerKernel, single(kspaceDst), option.performCenteredFFT, option.maxIterSPIRiT, option.stopThresSPIRiT, option.printInfo);
