
clear all
close all

cd C:\home\vessel_utilities\Ding_MR_Recon\MrRecon\compressSensing\SPIRiT_v0.1\Release_v0.1

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

sparseTransform = RedundantHarrDWT2D(1);
coeff = sparseTransform*im;
[m, totalNorm] = sparseTransform.coeffNorm(coeff);
totalNorm

norm(im(:))

im2 = sparseTransform'*coeff;
norm(im2(:)-im(:))

coeff = Matlab_PerformRedundantHarrDWT2D(im, 1, 1);
im2 = Matlab_PerformRedundantHarrDWT2D(coeff, 1, -1);
imtool(real([im(:,:,1) im2(:,:,1)]), []);
imtool(imag([im(:,:,1) im2(:,:,1)]), []);

% -----------------------------------------------------

sparseTransform2 = UndecimatedWavelet(1, 'db1', 'ppd');
coeff2 = sparseTransform2*im;
[m2, totalNorm2] = sparseTransform2.coeffNorm(coeff2);
totalNorm2

im2 = sparseTransform2'*coeff2;
norm(im2(:)-im(:))

% -----------------------------------------------------

imtool(coeff2{1,1}.dec{2}, []);
imtool(m2{1}, []);

coeff3 = cell(8,1);
for c=1:8
    coeff3{c,1}.dec = cell(4,1);
    coeff3{c,2}.dec = cell(4,1);
    
    for w=1:4
        coeff3{c,1}.dec{w} = real(coeff(:,:,w,c));
        coeff3{c,2}.dec{w} = imag(coeff(:,:,w,c));
    end
    
    coeff3{c,1}.sizeINI = [200 200];
    coeff3{c,2}.sizeINI = [200 200];
end
[m3, totalNorm3] = sparseTransform2.coeffNorm(coeff3);
totalNorm3

norm(m(:,:,1)-m3{1})
norm(m(:,:,2)-m3{2})
norm(m(:,:,3)-m3{3})

% compute the diff on wavelet coeff
for c=1:1
    for w=1:4
        myWav = coeff(:,:,w,c);
        mWav = coeff2{c,1}.dec{w}+i*coeff2{c,2}.dec{w};
        mWav2 = crop(mWav, [200 200]);
        norm(myWav(:))
        % norm(mWav(:))
        norm(mWav2(:))
    end
end

for w=1:3
    myNorm = m(:,:,1);
    mNorm = m2{w};
    mNorm2 = crop(mNorm, [200 200]);
    norm(myNorm(:))
    norm(mNorm2(:))
    imtool([myNorm mNorm2], [])
end

impad = zpad(im,256,256,8);
res= Matlab_PerformCudaDWT2D(single(impad), 1,1);
imtool(abs(res(:,:,1)), []);
imtool(abs(res(:,:,2)), []);
imtool(abs(res(:,:,3)), []);

Matlab_SaveAnalyze(single(real(im)), header, 'c:/real2.hdr');
Matlab_SaveAnalyze(single(imag(im)), header, 'c:/imag2.hdr');

header = CreateFtkHeaderInfo(impad, [1 1 1])
Matlab_SaveAnalyze(single(real(impad)), header, 'c:/real.hdr');
Matlab_SaveAnalyze(single(imag(impad)), header, 'c:/imag.hdr');

[r, header] = Matlab_LoadAnalyze('c:/realWav.hdr');
[iiii, header] = Matlab_LoadAnalyze('c:/imagWav.hdr');

pp = r + i*iiii;

im3 = ifft2c(pp);
plotComplexImage(im3, [1 1 1], 1204, 1024);

norm(im3(:)-im(:))
