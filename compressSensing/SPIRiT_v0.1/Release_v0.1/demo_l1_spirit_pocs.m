
clear all
close all

cd D:\vessel_utilities\Ding_MR_Recon\MrRecon\compressSensing\SPIRiT_v0.1\Release_v0.1

%load phantom.mat
load brain_8ch


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

% ------------------------------------------------------

% sparseTransform = RedundantHarrDWT2D(2);
% coeff = sparseTransform*im;
% [m, totalNorm] = sparseTransform.coeffNorm(coeff);
% totalNorm
% 
% im2 = sparseTransform'*coeff;
% norm(im2(:)-im(:))
% 
% coeff = Matlab_PerformRedundantHarrDWT2D(im, 1, 1);
% im2 = Matlab_PerformRedundantHarrDWT2D(coeff, 1, -1);
% imtool(real([im(:,:,1) im2(:,:,1)]), []);
% imtool(imag([im(:,:,1) im2(:,:,1)]), []);
% 
% % -----------------------------------------------------
% 
% sparseTransform2 = UndecimatedWavelet(1, 'db1', 'ppd');
% coeff2 = sparseTransform2*im;
% [m2, totalNorm2] = sparseTransform2.coeffNorm(coeff2);
% totalNorm2
% 
% im2 = sparseTransform2'*coeff2;
% norm(im2(:)-im(:))
% 
% % -----------------------------------------------------
% 
% imtool(coeff2{1,1}.dec{2}, []);
% imtool(m2{1}, []);
% 
% coeff3 = cell(8,1);
% for c=1:8
%     coeff3{c,1}.dec = cell(4,1);
%     coeff3{c,2}.dec = cell(4,1);
%     
%     for w=1:4
%         coeff3{c,1}.dec{w} = real(coeff(:,:,w,c));
%         coeff3{c,2}.dec{w} = imag(coeff(:,:,w,c));
%     end
%     
%     coeff3{c,1}.sizeINI = [200 200];
%     coeff3{c,2}.sizeINI = [200 200];
% end
% [m3, totalNorm3] = sparseTransform2.coeffNorm(coeff3);
% totalNorm3
% 
% norm(m(:,:,1)-m3{1})
% norm(m(:,:,2)-m3{2})
% norm(m(:,:,3)-m3{3})
% 
% % compute the diff on wavelet coeff
% for c=1:1
%     for w=1:4
%         myWav = coeff(:,:,w,c);
%         mWav = coeff2{c,1}.dec{w}+i*coeff2{c,2}.dec{w};
%         mWav2 = crop(mWav, [200 200]);
%         norm(myWav(:))
%         % norm(mWav(:))
%         norm(mWav2(:))
%     end
% end
% 
% for w=1:3
%     myNorm = m(:,:,1);
%     mNorm = m2{w};
%     mNorm2 = crop(mNorm, [200 200]);
%     norm(myNorm(:))
%     norm(mNorm2(:))
%     imtool([myNorm mNorm2], [])
% end

% impad = zpad(im,256,256,8);
% res= Matlab_PerformCudaDWT2D(single(impad), 1,1);
% imtool(abs(res(:,:,1)), []);
% imtool(abs(res(:,:,2)), []);
% imtool(abs(res(:,:,3)), []);
% 
% Matlab_SaveAnalyze(single(real(im)), header, 'c:/real2.hdr');
% Matlab_SaveAnalyze(single(imag(im)), header, 'c:/imag2.hdr');
% 
% header = CreateFtkHeaderInfo(impad, [1 1 1])
% Matlab_SaveAnalyze(single(real(impad)), header, 'c:/real.hdr');
% Matlab_SaveAnalyze(single(imag(impad)), header, 'c:/imag.hdr');
% 
% [r, header] = Matlab_LoadAnalyze('c:/realWav.hdr');
% [iiii, header] = Matlab_LoadAnalyze('c:/imagWav.hdr');
% 
% pp = r + i*iiii;
% 
% im3 = ifft2c(pp);
% plotComplexImage(im3, [1 1 1], 1204, 1024);
% 
% norm(im3(:)-im(:))
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

% tic
% [AtA,] = corrMatrix(kCalib,kSize);
% for n=1:coils
% 	kernel(:,:,:,n) = calibrate(AtA,kSize,coils,n,CalibTyk);
% end
% toc
% GOP = SPIRiT(kernel, 'fft',[fe,pe]);

tic 
kernel = spiritCalibration(kCalib, kSize, CalibTyk);
toc
GOP = SPIRiT(kernel, 'fft',[fe,pe]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                  Reconstruction                        %%%%%%%%%


disp('performing pocs reconstruction')
tic;
[res_pocs, prevXs] = pocsSPIRiT(DATA,GOP,nIter,DATA,wavWeight,1, 0.0022, 1024, 1024);
% [res, RESVEC] = cgSPIRiT_GPU(DATA,GOP, nIter, 1e-4, DATA);
toc

im_cgSPIRiT = ifft2c(res);
diffCGSPRiRIT = norm(im_cgSPIRiT(:)-im(:));

figure
hold on
plot(rAll, Diff, 'r-');
plot(rAll, Diff, 'bx');
plot(rAll(1):rAll(end), diffCGSPRiRIT, '-');
hold off

plotKSpace(res_pocs, [1 1 1], 1204, 1024);
plotKSpace(res, [1 1 1], 1204, 1024);

im_pocsspirit = ifft2c(res_pocs);
im_grappa = ifft2c(res_grappa);

im_pocsspirit_err = im_pocsspirit - im;
im_grappa_err = im_grappa - im;

N = numel(prevXs);
nRMSE = zeros(N,1);
pixelN = numel(im(:));
maxD = max(abs(im(:)));
minD = min(abs(im(:)));

for k=1:N
    imRecon = ifft2c(prevXs{k});
    diff = im(:)-imRecon(:);
    nRMSE(k) = sqrt(sumsqr(abs(im(:)-imRecon(:)))/pixelN)/(maxD-minD);
end

figure; 
hold on 
plot(nRMSE);
plot(nRMSE, 'x');
hold off

%% regularized CG

TVWeight = 0; 	% Weight for TV penalty
xfmWeight = 0.0025;	% Weight for Transform L1 penalty
Itnlim = 40;		% Number of iterations

params = initRegularizedCGParams();
params.TVWeight =TVWeight;     % TV penalty 
params.xfmWeight = xfmWeight;  % L1 wavelet penalty
params.Itnlim = Itnlim;

params.show = 1;

params.sparseTransform = UndecimatedWavelet(2, 'db1', 'ppd');

params

% disp('performing regularized CG reconstruction')
% tic;
% res_regCG = cgRegularizedSPIRiT(DATA, GOP, Itnlim, DATA, params);
% toc
% 
% im_regCG = ifft2c(res_regCG);

params.sparseTransform = UndecimatedWavelet(2, 'db1', 'ppd');
tic;
res_regCG2 = cgRegularizedSPIRiTWithLineSearch(DATA, GOP, Itnlim, DATA, params);
toc
im_regCG2 = ifft2c(res_regCG2);
plotComplexImage(im_regCG2, [1 1 1], 1204, 1024);


GOP1 = SPIRiT(kernel, 'fft',[fe,pe]);
params.sparseTransform = RedundantHarrDWT2D(2);
tic;
res_regCG2 = cgRegularizedSPIRiTWithLineSearch(DATA, GOP1, Itnlim, DATA, params);
toc
im_regCG2 = ifft2c(res_regCG2);
plotComplexImage(im_regCG2, [1 1 1], 1204, 1024);

GOP2 = SPIRiT(kernel, 'fft_AllGPU',[fe,pe]);
params.sparseTransform = RedundantHarrDWT2D(2);
tic;
res_regCG2 = cgRegularizedSPIRiTWithLineSearch_GPU(DATA, GOP2, Itnlim, DATA, params);
toc

[combinedImSen, sen] = adaptiveCoilCombination2DForSensitivity_PKMethod(im_regCG2, 32);   
plotComplexImageArrayMontage(sen, [2.44 2.44 10], 2048, 2048, 2, 1, 1/24);
plotComplexImage(combinedImSen, [1 1 1], 1204, 1024);

diffSPRiRIT = norm(im_regCG2(:)-im(:))

save brane_8ch_regGrappaTest

%% PRUNO recon
rAll = 10:10:180;
WdAll = [5 7 9 11 13];
resPURNOAll = zeros([size(DATA) numel(rAll) numel(WdAll)]);

for w=1:numel(WdAll)
    w
    for rr=1:numel(rAll)
        rr    
        r = rAll(rr);
        Wd = WdAll(w);
        
        pru = PRUNO(Wd, r);

        thresEigenValue = 0.01
        tic
        kernel = pru.prunoCalibration(kCalib,[Wd Wd],-1, r, 0);
        toc

        tic
        kernelIm = pru.convertKernelToImage(kernel, [fe,pe], 0, [1 1 1], 1024, 1024);
        toc

        tic
        kernelComposite = pru.computeCompositeNullKernel(kernelIm, 0, [1 1 1], 1024, 1024);
        toc

        pru.kernel = kernel;
        pru.kernelIm = kernelIm;
        pru.kernelComposite = kernelComposite;

        [resPURNO, RESVEC] = cgPRUNO(DATA, pru, 300, 1e-6, DATA);
    %     im_PURNO = ifft2c(resPURNO);
    %     plotComplexImage(im_PURNO, [1 1 1], 1204, 1024);

        [resPURNO, RESVEC] = cgPRUNO_NotLSQR(DATA, pru, 300, 1e-9, DATA);
        resPURNOAll(:,:,:,rr,w) = resPURNO;
    end
end

Diff = zeros([numel(rAll) numel(WdAll)]);

for w=1:numel(WdAll)
    for rr=1:numel(rAll)
         im_PURNO = ifft2c(resPURNOAll(:,:,:,rr,w));
    %     plotComplexImage(im_PURNO, [1 1 1], 1204, 1024);

        Diff(rr,w) = norm(im_PURNO(:)-im(:));
    end
end

figure
hold on
plot(rAll, Diff(:,1), 'r-');
plot(rAll, Diff(:,2), 'b-');
plot(rAll, Diff(:,3), 'y-');
plot(rAll, Diff(:,4), 'c-');
plot(rAll, Diff(:,5), 'g-');
plot(rAll(1):1:rAll(end), diffCGSPRiRIT*ones(size(rAll(1):1:rAll(end))), 'k');
hold off
legend('w=5','w=7','w=9','w=11','w=13','SPIRiT');
ylabel('norm(groundtruth-recon)');
xlabel('number of kernels r');

% plotComplexImageArrayMontage(ifft2c(resPURNOAll), [1 1 1], 1024, 1024, 1, 1, 0);

plotComplexImage(ifft2c(resPURNOAll(:,:,:,16,1)), [1 1 1], 1204, 1024);

save PURNO_Change_r_w

% ------------------------------------------------------------------
params.xfmWeight = 0.001;
params.grappaPriorWeight = 0.1;
tic;
res_regCGSen = cgRegularizedSPIRiTWithLineSearchAndGrappaPrior(DATA, GOP, sen, res_regCG2, Itnlim, DATA, params);
toc
im_regCGSen = ifft2c(res_regCGSen);
plotComplexImage(im_regCGSen, [1 1 1], 1204, 1024);

% ------------------------------------------------------------------

params.xfmWeight = 0.01
tic
res_regGrappa = cgRegularizedGRAPPAWithLineSearch(DATA, res_regCG2, sen, Itnlim, res_regCG2, params);
toc
im_regGrappa = ifft2c(res_regGrappa);
plotComplexImage(im_regGrappa, [1 1 1], 1204, 1024);

% NESTA optimization
% do NESTA recon

% params.sparseTransform = Wavelet('Haar',2,3);
% res_nesta = nestaSPIRiT(DATA, GOP, Itnlim, DATA, params);
% im_nesta = ifft2c(res_nesta);

params.sparseTransform = Wavelet('Haar',2,3);
params.centre = 1024;
params.width = 1024;
res_fista = fistaSPIRiT(DATAcomp, GOP, Itnlim, DATAcomp, params);
im_fista = ifft2c(res_fista);
plotComplexImage(im_fista, [1 1 1], 1204, 1024);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%		Display                                  %%%%%%%%%%%

plotComplexImage(im, [1 1 1], 1204, 1024);
plotComplexImage(im_dc, [1 1 1], 1204, 1024);
plotComplexImage(im_grappa, [1 1 1], 1204, 1024);
plotComplexImage(im_pocsspirit, [1 1 1], 1204, 1024);
plotComplexImage(im_regCG, [1 1 1], 1204, 1024);

im_regCG_sqr = sos(im_regCG);
im_pocsspirit_sqr = sos(im_pocsspirit);
im_grappa_sqr = sos(im_grappa);
im_dc_sqr = sos(im_dc);
im_sqr = sos(im);
im_pocsspirit_err_sqr = sos(im_pocsspirit_err);
im_grappa_err_sqr = sos(im_grappa_err);


figure, imshow(cat(1,cat(2,im_sqr,im_dc_sqr, zeros(size(im_sqr))),cat(2,im_pocsspirit_sqr,im_grappa_sqr, im_regCG_sqr)),[],'InitialMagnification',150);
title ('top-left:Full		top-right:zero-fill w/dc	bottom-left:SPIRiT 	 	bottom-right:GRAPPA');
figure, imshow(cat(2,im_pocsspirit_err_sqr,im_grappa_err_sqr),[],'InitialMagnification',150);
title ('Difference images: SPIR-iT               GRAPPA');
