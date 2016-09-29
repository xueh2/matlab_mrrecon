
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
plotKSpace(DATA);

R = 2;
Nfe = size(DATA, 1)
Npe = size(DATA, 2)
CHA = size(DATA, 3)

undersampledKspace = zeros(size(DATA));
undersampledKspace(:, 1:R:Npe, :) = DATA(:, 1:R:Npe, :);
plotKSpace(undersampledKspace);

coilMapMethod = 'Souheil';
dstChaThres = 0.001;
performSENSE = 0;
performGRAPPA = 1;
performVICPAAS = 0;

option = CreateMrReconOption(R, R, coilMapMethod, dstChaThres, performSENSE, performGRAPPA, performVICPAAS, Nfe, Nfe/2, Npe, Npe/2, Npe/2);

option.acqAcsLines = 1:Npe;
option.acqPELines = 1:Npe;

% Coef = GRAPPA_SrcDstChannels_Kernel_2D(DATA, DATA, option, option.thresReg);
% kernelIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, option, Nfe, Npe, R);

[refFull, E_0, V, dstCha] = coilReduction(DATA, 7, -1);
dst = applyEigenVector(DATA, V);
dst = dst(:,:,end-dstCha+1:end,:);

complexImageSen = ifft2c(dst);
sensitivityMap = CoilSensitivity_Souheil_Walsh(complexImageSen, [7 7], 0);

Coef = GRAPPA_SrcDstChannels_Kernel_2D_2(DATA, dst, option, option.thresReg);
kernelIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients_2(Coef, option, Nfe, Npe, R);
Coef3D = reshape(Coef, [5 4 1 8 8 2 1]);

unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff(sensitivityMap, kernelIm);
unWarppedCombinedIm = HTGRAPPA_SrcDstChannels_ImagedomainReconWithUnmixCoeff(undersampledKspace, unmixCoeff, []);
imagescn(abs(unWarppedCombinedIm));

option.KernelSize(3) = 1;
option.KernelPatternPAR = [0];
option.OutPatternPAR = [0];

Coef3D = GRAPPA_SrcDstChannels_Kernel_3D(DATA, dst, option, option.thresReg);
kernelIm3D = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients_3D(Coef3D, option, Nfe, Npe, 1, R, 1);

unmixCoeff3D = HTGRAPPA_SrcDstChannels_UnmixingCoeff_3D(sensitivityMap, kernelIm3D);
figure; imagescn(abs(unmixCoeff3D));

unwrappedIm3D = HTGRAPPA_SrcDstChannels_ImagedomainReconWithUnmixCoeff_3D(undersampledKspace, unmixCoeff3D);
imagescn(abs(unwrappedIm3D));

size(kernelIm3D)
kernelIm = squeeze(kernelIm3D);

aliasedIm = ifft2c(undersampledKspace);
aliasedIm = repmat(aliasedIm, [1 1 1 dstCha]);
unWarppedIm = sum(aliasedIm.*kernelIm, 3);
unWarppedIm = squeeze(unWarppedIm);
imagescn(SoS_Image(unWarppedIm));

