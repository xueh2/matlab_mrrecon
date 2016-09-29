
clear all
close all

cd F:\data\OSU_RO1_Grant\CompressSensing\FreeBreathingCine_20110531

% load('Free Breathing Cine.mat')
% fullKSpace = performDownSampleFE(Cine_12ch_free);

load fullKSpace

centre = 250
width = 500

plotKSpaceArray(fullKSpace, [1 1 1], centre, width, 1/24);

DATA = fullKSpace(:,:,:,20);
DATA = double(DATA);
plotKSpace(DATA, [1 1 1], centre, width, 1, 1/24);

Im = ifft2DForVolume(DATA);
WT = UndecimatedWavelet(2, 'db1', 'ppd');
res = WT*Im;
D2 = WT'*res;
pp = WT.softThresh(res, 0.0015);

norm(D2(:)-Im(:))
plotComplexImage(Im, [1 1 1], centre, width, 1, 1/24);
plotComplexImage(D2, [1 1 1], centre, width, 1, 1/24);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Reconstruction Parameters %%%%%%%%%%%%%%%%%%%
%						     %	 			
kSize = [7 7];  % SPIRiT kernel size
nIter = 20; % number of iteration; phantom requires twice as much as the brain.
CalibTyk = 1e-2;  % Tykhonov regularization in the calibration
wavWeight = 0.0015;  % Wavelet soft-thresholding regularization in the reconstruction (SPIRiT only)
skipGRAPPA = 1; % skip GRAPPA recon.
showFlag = 1;

im = ifft2c(DATA);
plotComplexImage(im, [1 1 1], centre, width, 1, 1/24);

% generate variable density random sampling
N = size(DATA); 		% image Size
DN = size(DATA); 	% data Size
pctg = [0.25];  	% undersampling factor
P = 5;			% Variable density polymonial degree
Itnlim = 8;		% Number of iterations

pdf = genPDF(DN,P,pctg , 2 ,0.1,0);	% generates the sampling PDF
k = genSampling(pdf,10,60);		% generates a sampling pattern

% k = zeros(N(1), N(2));
% k(:, 1:4:end) = 1;

range = ceil(N(1)*0.1):floor(N(1)*0.9);
k(range, N(2)/2-10:N(2)/2+10) = 1;
mask = k;
imtool(mask, []);

[CalibSize, dcomp] = getCalibSize(mask);  % get size of calibration area from mask
pe = size(DATA,2); fe = size(DATA,1); coils = size(DATA,3); % get sizes
DATA = DATA.*repmat(mask,[1,1,coils]); % multiply with sampling matrix
plotKSpace(DATA, [1 1 1], centre, width, 1, 1/24);

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
plotComplexImage(im_dc, [1 1 1], centre, width, 1, 1/24);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GRAPPA                                   %%%%%%%%%%%
% if skipGRAPPA
%     disp('skipping grappa, replacing with zero filling');
%     res_grappa = DATA;
% else
%     disp('performing traditional GRAPPA reconstruction');
%     kCalib = crop(DATA,[CalibSize,coils]);
%     res_grappa = GRAPPA(DATA,kCalib,kSize,CalibTyk);
% end
% plotKSpace(res_grappa, [1 1 1], centre, width, 1, 1/24);
res_grappa = DATA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%		 Perform Calibration                       %%%%%%%%%%
disp('performing calibration for SPIRiT')
kCalib = crop(DATA,[CalibSize,coils]);
kernel = zeros([kSize,coils,coils]);

[AtA,] = corrMatrix(kCalib,kSize);
for n=1:coils
	kernel(:,:,:,n) = calibrate(AtA,kSize,coils,n,CalibTyk);
end
GOP = SPIRiT(kernel, 'fft',[fe,pe]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                  Reconstruction                        %%%%%%%%%

disp('performing pocs reconstruction')
tic;
[res_pocs, prevXs] = pocsSPIRiT(DATA,GOP,nIter,DATA,wavWeight,showFlag);
[res_pocs2, RESVEC] = cgSPIRiT(DATA,GOP, nIter, 1e-5, DATA);
toc
% [res_pocs] = pocsSPIRiT(DATA,GOP,nIter,DATA,0,1);

plotKSpace(res_pocs, [1 1 1], centre, width, 1, 1/24);

im_pocsspirit = ifft2c(res_pocs);
im_grappa = ifft2c(res_grappa);

im_pocsspirit_err = im_pocsspirit - im;
im_grappa_err = im_grappa - im;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%		Display                                  %%%%%%%%%%%
im_pocsspirit_sqr = sos(im_pocsspirit);
im_grappa_sqr = sos(im_grappa);
im_dc_sqr = sos(im_dc);
im_sqr = sos(im);
im_pocsspirit_err_sqr = sos(im_pocsspirit_err);
im_grappa_err_sqr = sos(im_grappa_err);


figure, imshow(cat(1,cat(2,im_sqr,im_dc_sqr),cat(2,im_pocsspirit_sqr,im_grappa_sqr)),[],'InitialMagnification',150);
title ('top-left:Full		top-right:zero-fill w/dc	bottom-left:SPIRiT 	 	bottom-right:GRAPPA');
figure, imshow(cat(2,im_pocsspirit_err_sqr,im_grappa_err_sqr),[],'InitialMagnification',150);
title ('Difference images: SPIR-iT               GRAPPA');
