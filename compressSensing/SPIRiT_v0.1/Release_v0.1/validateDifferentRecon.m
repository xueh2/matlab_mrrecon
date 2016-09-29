
clear all
close all

cd D:\vessel_utilities\mr\MrRecon\compressSensing\SPIRiT_v0.1\Release_v0.1

%load phantom.mat
load brain_8ch

DATAOri = DATA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Reconstruction Parameters %%%%%%%%%%%%%%%%%%%
%						     %	 			
kSize = [5,5];  % SPIRiT kernel size
nIter = 20; % number of iteration; phantom requires twice as much as the brain.
mask_type = 'random4'; % options are: 'unif','random4','random3'
CalibTyk = 0.01;  % Tykhonov regularization in the calibration
wavWeight = 0.0015;  % Wavelet soft-thresholding regularization in the reconstruction (SPIRiT only)
addNoiseSTD = 0.0; % add noise. Use only for phantom. Brain has enough noise!
skipGRAPPA = 1; % skip GRAPPA recon.

DATA = DATA+randn(size(DATA))*addNoiseSTD + i*randn(size(DATA))*addNoiseSTD;
im = ifft2c(DATA);

switch mask_type
    case 'unif'
           mask = mask_unif;
           disp('Using uniform undersampling.')
           disp('Change mask_type in the code for other sampling patterns.');
           disp(' ');
           disp(' ');
           
    case 'random3'
            mask = mask_randm_x3;
            if skipGRAPPA==0
                disp('Warning:GRAPPA recon is slow for random sampling.')
                disp('change to skipGRAPPA=1 in the code to skipp.')
                disp(' ');
                disp(' ');
                
                gSkip=1;
            end
            
            
    case 'random4'
            mask = mask_randm_x4;
            if skipGRAPPA==0
                disp('Warning:GRAPPA recon is slow for random sampling.')
                disp('change to skipGRAPPA=1 in the code to skipp.')
                disp(' ')
                disp(' ')
                
                gSkip=1;
            end
    otherwise
        mask = mask_unif'
end

[CalibSize, dcomp] = getCalibSize(mask);  % get size of calibration area from mask
pe = size(DATA,2); fe = size(DATA,1); coils = size(DATA,3); % get sizes
DATA = DATA.*repmat(mask,[1,1,coils]); % multiply with sampling matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale the data such that the zero-filled density compensated      %%%%%%%%%
% k-space norm is 1. This is useful in order to use similar         %%%%%%%%%
% regularization penalty values for different problems.             %%%%%%%%%

DATAcomp = DATA.*repmat(dcomp,[1,1,coils]);
scale_fctr = norm(DATAcomp(:))/sqrt(coils)/20;
DATA = DATA/scale_fctr;
DATAcomp = DATAcomp/scale_fctr;

im_dc = ifft2c(DATAcomp);
im = im/scale_fctr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GRAPPA                                   %%%%%%%%%%%
if skipGRAPPA
    disp('skipping grappa, replacing with zero filling');
    res_grappa = DATA;
else
    disp('performing traditional GRAPPA reconstruction');
    kCalib = crop(DATA,[CalibSize,coils]);
    res_grappa = GRAPPA(DATA,kCalib,kSize,CalibTyk);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%		 Perform Calibration                       %%%%%%%%%%
disp('performing calibration for SPIRiT')
kCalib = crop(DATA,[CalibSize,coils]);
kernel = zeros([kSize,coils,coils]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                  Reconstruction                        %%%%%%%%%
loadGT = 0;
if ( isFileExist('validationReconAll.mat') )
    load('validationReconAll.mat');
    res_pocs_V = res_pocs;
    res_regCG_UndecimatedWavelet_V = res_regCG_UndecimatedWavelet;
    res_regCG_RedundantHarrDWT2D_V = res_regCG_RedundantHarrDWT2D;
    res_regCG_RedundantHarrDWT2D_GPU_V  = res_regCG_RedundantHarrDWT2D_GPU;
    res_regCG_Harr_NoNullSpace_V  = res_regCG_Harr_NoNullSpace;
    res_regCG_Harr_NoNullSpace_GPU_V  = res_regCG_Harr_NoNullSpace_GPU;
    loadGT = 1;
end


tic 
kernel = spiritCalibration(kCalib, kSize, CalibTyk);
toc
GOP = SPIRiT(kernel, 'fft',[fe,pe]);

disp('performing pocs reconstruction')
tic;
[res_pocs, prevXs] = pocsSPIRiT(DATA,GOP,nIter,DATA,wavWeight,0, 0.0022, 1024, 1024);
toc
plotComplexImage(ifft2c(res_pocs), [1 1 1], 1204, 1024);

%% regularized CG

TVWeight = 0; 	% Weight for TV penalty
xfmWeight = 0.0025;	% Weight for Transform L1 penalty
Itnlim = 40;		% Number of iterations

params = initRegularizedCGParams();
params.TVWeight =TVWeight;     % TV penalty 
params.xfmWeight = xfmWeight;  % L1 wavelet penalty
params.Itnlim = Itnlim;
params.show = 0;

params.sparseTransform = UndecimatedWavelet(2, 'db1', 'ppd');
tic;
res_regCG_UndecimatedWavelet = cgRegularizedSPIRiTWithLineSearch(DATA, GOP, Itnlim, DATA, params);
toc

GOP1 = SPIRiT(kernel, 'fft',[fe,pe]);
params.sparseTransform = RedundantHarrDWT2D(2);
tic;
[kspaceLinear, RESVEC] = cgSPIRiT(DATA,GOP1,Itnlim,1e-4,DATA);
res_regCG_RedundantHarrDWT2D = cgRegularizedSPIRiTWithLineSearch(DATA, GOP1, Itnlim, kspaceLinear, params);
toc
plotComplexImage(ifft2c(kspaceLinear), [1 1 1], 1204, 1024);
plotComplexImage(ifft2c(res_regCG_RedundantHarrDWT2D), [1 1 1], 1204, 1024);

GOP2 = SPIRiT(kernel, 'fft_AllGPU',[fe,pe]);
params.sparseTransform = RedundantHarrDWT2D(2);
tic;
[kspaceLinear, RESVEC] = cgSPIRiT_GPU(DATA,GOP2,Itnlim,1e-4,DATA);
res_regCG_RedundantHarrDWT2D_GPU = cgRegularizedSPIRiTWithLineSearch_GPU(DATA, GOP2, Itnlim, kspaceLinear, params);
toc
plotComplexImage(ifft2c(kspaceLinear), [1 1 1], 1204, 1024);
plotComplexImage(ifft2c(res_regCG_RedundantHarrDWT2D_GPU), [1 1 1], 1204, 1024);

[kspaceLinear, iterFinal] = cgSPIRiT_NotLSQR_GPU(DATA,GOP2,Itnlim,1e-9,zeros(size(DATA)));
plotComplexImage(ifft2c(kspaceLinear), [1 1 1], 1204, 1024);

params.dataWeight = 0.5;
params.sparseTransform = RedundantHarrDWT2D(2);
tic
[kspaceLinear, RESVEC] = cgSPIRiT(DATA,GOP1,Itnlim,1e-4,DATA);
res_regCG_Harr_NoNullSpace = cgRegularizedSPIRiTWithLineSearch_WithoutNullSpace(DATA, GOP1, Itnlim, kspaceLinear, params);
toc

params.dataWeight = 0.5;
params.sparseTransform = RedundantHarrDWT2D(2);
tic
[kspaceLinear, RESVEC] = cgSPIRiT_GPU(DATA,GOP2,Itnlim,1e-4,DATA);
res_regCG_Harr_NoNullSpace_GPU = cgRegularizedSPIRiTWithLineSearch_WithoutNullSpace_GPU(DATA, GOP2, Itnlim, kspaceLinear, params);
toc

GOPRep = SPIRiT_KernelIm(GOP2.KERNEL, 'fft_AllGPU');
[kspaceLinear2, RESVEC] = cgSPIRiT_GPU(DATA,GOPRep,Itnlim,1e-4,DATA);
plotComplexImage(ifft2c(kspaceLinear2), [1 1 1], 1204, 1024);

% vicpaas
[DATAOriDst, E_0, V, dstCha] = coilReduction(DATAOri, -1, -1);
kspaceDst = applyEigenVector(DATA, V);

figure; imagescn(abs(ifft2c(DATAOri)));
figure; imagescn(abs(ifft2c(DATAOriDst)));

s = size(DATAOri);
p = reshape(DATAOri, [s(1)*s(2) s(3)]);
p2 = p*V;
p3 = reshape(p2, s);
figure; imagescn(abs(ifft2c(p3)));

p4 = p2*V';
p4 = reshape(p4, s);
figure; imagescn(abs(ifft2c(p4)));

imSize = [size(DATAOri,1) size(DATAOri,2)];

optionVSpirit.KernelSize = [7 7];
optionVSpirit.OutPattern = [1 1];
kernelS2D = spiritCalibration_SrcDst_MultipleOutput(DATAOri, DATAOriDst, optionVSpirit, CalibTyk);
kernelD2S = spiritCalibration_SrcDst_MultipleOutput(DATAOriDst, DATAOri, optionVSpirit, CalibTyk);    

KERNELS2D = convertSPIRiTKSpaceKernel2ImageSpace(kernelS2D, [fe pe]);
KERNELD2S = convertSPIRiTKSpaceKernel2ImageSpace(kernelD2S, [fe pe]);

GP = spiritImageDomainKernelCombined(KERNELS2D, KERNELD2S);

kfe = size(kernelS2D, 1);
kpe = size(kernelS2D, 2);

hkfe = ceil(kfe/2);
hkpe = ceil(kpe/2);
            
cha = size(DATAOriDst, 3);
identityKernel = zeros(kfe,kpe,cha,cha);
for d=1:cha
    identityKernel(hkfe, hkpe, d, d) = 1.0;
end
identityKernelIm = convertSPIRiTKSpaceKernel2ImageSpace(identityKernel, [fe pe]);

GP_I = GP - identityKernelIm;

VIC = VICPAAS_GP_I(GP_I, 'fft_AllGPU');
[res, RESVEC] = cgSPIRiT_GPU(kspaceDst, VIC, 100,1e-4, kspaceDst);
plotComplexImage(ifft2c(res), [1 1 1], 1204, 1024);

gpuGOPs = cell(3,1);
gpuGOPs{1} = VIC;
gpuGOPs{2} = VIC;
gpuGOPs{3} = VIC;

data2 = zeros([size(kspaceDst) 3]);
data2(:,:,:,1) = kspaceDst;
data2(:,:,:,2) = kspaceDst;
data2(:,:,:,3) = kspaceDst;

GOP2DT = SPIRiT2DPlusT(gpuGOPs);

dd = GOP2DT'*data2;
ddA = VIC'*kspaceDst;

ddDiff = dd(:,:,:,1)-ddA;
norm(ddDiff(:))
norm(gather(ddDiff(:)))

[res, RESVEC] = cgSPIRiT2DPlusT_GPU(data2, GOP2DT, Itnlim,1e-4,data2);
plotComplexImage(ifft2c(res(:,:,:,1)), [1 1 1], 1204, 1024);
plotComplexImage(ifft2c(res(:,:,:,2)), [1 1 1], 1204, 1024);

cgRegularizedSPIRiTWithLineSearch_2DPlusT_CoilSensitivity_GPU

VIC = VICPAAS(kernelS2D, kernelD2S, 'fft_AllGPU',[fe,pe], V);
params.sparseTransform = RedundantHarrDWT2D(2);
tic;
[kspaceLinear, RESVEC] = cgSPIRiT_GPU(DATA,VIC,100,1e-4,DATA);
res_regCG_RedundantHarrDWT2D_GPU = cgRegularizedSPIRiTWithLineSearch_GPU(DATA, GOP2, Itnlim, kspaceLinear, params);
toc
plotComplexImage(ifft2c(kspaceLinear), [1 1 1], 1204, 1024);
plotComplexImage(ifft2c(res_regCG_RedundantHarrDWT2D_GPU), [1 1 1], 1204, 1024);

GP_I = SPIRiTSrcDstCombined(kernelS2D,kernelD2S, 'fft_AllGPU', [fe pe]);
[res, RESVEC] = cgSPIRiT_SrcDstCombined_GPU(kspaceDst, GP_I, 100,1e-4, kspaceDst);

refSrc = DATAOri;
headerSrc = CreateFtkHeaderInfo(refSrc, [1 1 1 1]);
norm(refSrc(:))

refDst = DATAOriDst;
headerDst = CreateFtkHeaderInfo(refDst, [1 1 1 1]);
norm(refDst(:))

Matlab_SaveAnalyze(real(single(refSrc)), headerSrc, [ 'refSrc_real.hdr']);
Matlab_SaveAnalyze(imag(single(refSrc)), headerSrc, [ 'refSrc_imag.hdr']);

Matlab_SaveAnalyze(real(single(refDst)), headerDst, [ 'refDst_real.hdr']);
Matlab_SaveAnalyze(imag(single(refDst)), headerDst, [ 'refDst_imag.hdr']);

refDst = refSrc;

[kernelS2D, kernelD2S, kernelImS2D, conjKernelImS2D, kernelImD2S, conjKernelImD2S, GP_Im, identityKernelIm, GP_I_Im, conjGP_I_Im, identityKernel] = Matlab_PerformSPIRiTKernelCalibrationSrcDst2D(single(refSrc), headerSrc, single(refDst), headerDst, ... 
                    optionVSpirit.KernelSize, optionVSpirit.OutPattern, [fe pe], CalibTyk);
       
VIC = VICPAAS_GP_I(GP_I_Im, 'fft_AllGPU');

[res, RESVEC] = cgSPIRiT_GPU(kspaceDst, VIC, 100,1e-4, kspaceDst);


GOP2 = SPIRiT_KernelIm(kernelImS2D, 'fft_AllGPU');
[kspaceLinear, RESVEC] = cgSPIRiT_GPU(kspaceFrame, GOP2, 200, 1e-6, kspaceFrame);

kr = analyze75read('D:\software\Gadgetron\20130114\gadgetron\toolboxes\gtplus\ut\result\kerS2D_REAL.hdr');
ki = analyze75read('D:\software\Gadgetron\20130114\gadgetron\toolboxes\gtplus\ut\result\kerS2D_IMAG.hdr');
ker = kr + i*ki;

cd D:\software\Gadgetron\20130114\gadgetron\toolboxes\gtplus\ut\result
ker2 = permute(ker, [2 1 3 4]);
ker2 = single(ker2);
header = CreateFtkHeaderInfo(ker2, [1 1 1 1]);
Matlab_SaveAnalyze(real(ker2), header, 'kerS2D_REAL_2.hdr');
Matlab_SaveAnalyze(imag(ker2), header, 'kerS2D_IMAG_2.hdr');

figure; imagescn(abs(kernelS2D), [], [], [], 4);
figure; imagescn(abs(ker), [], [], [], 4);
figure; imagescn(abs(cat(5, ker, kernelS2D)), [], [], [], 5);


kImS2D = readGTPlusExportData('D:\software\Gadgetron\20130114\gadgetron\toolboxes\gtplus\ut\result\kerImS2D');
figure; imagescn(abs(kImS2D), [], [], [], 4);

kImS2D2 = readGTPlusExportData('D:\software\Gadgetron\20130114\gadgetron\toolboxes\gtplus\ut\result\kernelS2D_new_Im.hdr');
figure; imagescn(abs(kImS2D2), [], [], [], 4);

cd D:\software\Gadgetron\20130114\gadgetron\toolboxes\gtplus\ut\result

[data1, data2] = compareGTPlusExportData('kernelS2D_new_Im.hdr', 'kerImS2D');

[data1, data2] = compareGTPlusExportData('concatenatedKim.hdr', 'adjForwardkerIm');

[data1, data2] = compareGTPlusExportData('concatenatedKim.hdr', 'adjForwardkerIm');


%% 2D + T spirit
GOPs = cell(3,1);
GOPs{1} = GOP1;
GOPs{2} = GOP1;
GOPs{3} = GOP1;

gpuGOPs = cell(3,1);
gpuGOPs{1} = GOP2;
gpuGOPs{2} = GOP2;
gpuGOPs{3} = GOP2;

data2 = zeros([size(DATA) 3]);
data2(:,:,:,1) = DATA;
data2(:,:,:,2) = DATA;
data2(:,:,:,3) = DATA;

GOP2DT = SPIRiT2DPlusT(gpuGOPs);

dd = GOP2DT'*DATA;

ddA = GOP2'*DATA;

ddDiff = dd(:,:,:,1)-ddA;
norm(ddDiff(:))
norm(gather(ddDiff(:)))

GOPs = cell(1,1);
GOPs{1} = GOP1;
GOP2DT = SPIRiT2DPlusT(GOPs);
[res, RESVEC] = cgSPIRiT2DPlusT(DATA, GOP2DT, Itnlim,1e-4,DATA);
plotComplexImage(ifft2c(res(:,:,:,1)), [1 1 1], 1204, 1024);

gpuGOPs = cell(1,1);
gpuGOPs{1} = GOP2;

gpuGOP2DT = SPIRiT2DPlusT(gpuGOPs);
[res, RESVEC] = cgSPIRiT2DPlusT_GPU(data2, gpuGOP2DT, Itnlim,1e-4,data2);
plotComplexImage(ifft2c(res(:,:,:,1)), [1 1 1], 1204, 1024);
plotComplexImage(ifft2c(res(:,:,:,2)), [1 1 1], 1204, 1024);
plotComplexImage(ifft2c(res(:,:,:,3)), [1 1 1], 1204, 1024);

%% -----------------------------------------
plotComplexImage(ifft2c(res_pocs), [1 1 1], 1204, 1024);
plotComplexImage(ifft2c(res_regCG_UndecimatedWavelet), [1 1 1], 1204, 1024);
plotComplexImage(ifft2c(res_regCG_RedundantHarrDWT2D), [1 1 1], 1204, 1024);
plotComplexImage(ifft2c(res_regCG_RedundantHarrDWT2D_GPU), [1 1 1], 1204, 1024);
plotComplexImage(ifft2c(res_regCG_Harr_NoNullSpace_GPU), [1 1 1], 1204, 1024);

if ( loadGT )
    if ( norm(res_pocs_V(:)-res_pocs(:)) > 0 ) 
        error('error ... ')
    end
    
    diff = norm(res_regCG_UndecimatedWavelet_V(:)-res_regCG_UndecimatedWavelet(:))
    diff = diff + norm(res_regCG_RedundantHarrDWT2D_V(:)-res_regCG_RedundantHarrDWT2D(:))
    diff = diff + norm(res_regCG_RedundantHarrDWT2D_GPU_V(:)-res_regCG_RedundantHarrDWT2D_GPU(:))
    diff = diff + norm(res_regCG_Harr_NoNullSpace_V(:)-res_regCG_Harr_NoNullSpace(:))
    diff = diff + norm(res_regCG_Harr_NoNullSpace_GPU_V(:)-res_regCG_Harr_NoNullSpace_GPU(:))

    if ( diff > 0 ) 
        error('error ... ')
    end
        
    % save validationReconAll res_pocs res_regCG_UndecimatedWavelet res_regCG_RedundantHarrDWT2D res_regCG_RedundantHarrDWT2D_GPU res_regCG_Harr_NoNullSpace res_regCG_Harr_NoNullSpace_GPU

end

%% ------------------------
% validate 3D recon
load meas_MID91_t1_fl3d_tra_dynaVIEWS_p2p2_FID2026_ori
rx_dwelltime_data = 2900;

load meas_MID5799_T1_VIBE_FS_separate_2x2_FID50187_ori
rx_dwelltime_data = 3500;

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
thresReg = 0.01;
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
