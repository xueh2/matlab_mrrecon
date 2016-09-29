function [res, RESVEC] = cgPRUNO(y, PRU, nIter, lambda, x0)
%
%
%  res = cgPRUNO(y, PRU, nIter, lambda, x0)
%  
%  Implementation of the Cartesian, conjugated gradient PRUNO reconstruction
%  Solve PURNO equation: Dc*(N'N)*Dc'*d = -Dc*(N'N)*D'*d
%
%  Input:
%		y -	Undersampled k-space data. Make sure that empty entries are zero
%			or else they will not be filled.
%		PRU -	the PURNO operator obtained by calibration, implementing N'N
%		nIter -	Maximum number of iterations
%		lambda-	Tykhonov regularization parameter
%       x0    - Initil value
%
% Outputs:
%		res - Full k-space matrix
%
%
% Example:
%
%   [x,y] = meshgrid(linspace(0,1,128));
%   % Generate fake Sensitivity maps
%   sMaps = cat(3,x.^2,1-x.^2,y.^2,1-y.^2);
%   % generate 4 coil phantom
%   imgs = repmat(phantom(128),[1,1,4]).*sMaps;
%   DATA = fft2c(imgs);
%   % crop 20x20 window from the center of k-space for calibration
%   kCalib = crop(DATA,[20,20,4]);
%
%   %calibrate a kernel
%   kSize = [5,5];
%   coils = 4;
%   Wd = 8;
%   r = 50;
%   pru = PRUNO(Wd, r);
%   thresEigenValue = 0.01
%   kernel = pru.prunoCalibration(kCalib,kSize,thresEigenValue, r, 1);
% 
%   kernelIm = pru.convertKernelToImage(kernel, [fe,pe], 0, [1 1 1], 1024, 1024);
% 
%   kernelComposite = pru.computeCompositeNullKernel(kernelIm, 0, [1 1 1], 1024, 1024);
%   
%   res = PRUNO(Wd, r); 
%   res.kernel = kernel;
%   res.kernelIm = kernelIm;
%   res.kernelComposite = kernelComposite;
% 
%   % undersample by a factor of 2
%   DATA(1:2:end,2:2:end,:) = 0;
%   DATA(2:2:end,1:2:end,:) = 0;
%   
%   %reconstruct:
%   [res] = cgPRUNO(DATA, res, 20, 1e-5, DATA);
%   figure, imshow(cat(2,sos(imgs), 2*sos(ifft2c(DATA)), sos(ifft2c(res))),[]);
%   title('full,  zero-fill,   result')
%   
%   Hui Xue, July 2011
%

tCGPRUNO=tic;

if nargin < 4
	lambda = 0;
end

if nargin < 5
     x0 = y;
end

[sx,sy,nCoils] = size(y);

idx_acq = find(abs(y)>0);
idx_nacq = find(abs(y)==0);
N = length(idx_nacq(:));

% compute -Im(N'N)Ia*d 
% y is Ia*d, a kspace containing only acquired points
yy = PRU*y; % (N'N)*y 
yy(idx_acq) = 0; % apply Im, only nonacquired points are kept
yy = -yy(idx_nacq);

[tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,1e-6,nIter, speye(N,N),speye(N,N),x0(idx_nacq),PRU,sx,sy,nCoils,idx_acq,idx_nacq,lambda);
% [tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,1e-6,nIter, [],[],x0(idx_nacq),GOP,sx,sy,nCoils,idx_nacq, lambda);
FLAG
RELRES
ITER

res = y;
res(idx_nacq) = tmpres;

disp(['cgPRUNO : ' num2str(toc(tCGPRUNO))]);

function [res,tflag] = aprod(x,PRU,sx,sy,nCoils,idx_acq,idx_nacq,lambda,tflag)

	if strcmp(tflag,'transp');
		tmpx = zeros(sx,sy,nCoils);
		tmpx(idx_nacq) = x; % Im*d, get the nonacquired points, other points are zero
		res = PRU*tmpx;
        res(idx_acq) = 0; % apply another Im
        res = res(idx_nacq);	
    else
		tmpx = zeros(sx,sy,nCoils);
		tmpx(idx_nacq) = x; % Im*d, get the nonacquired points, other points are zero
		res = PRU*tmpx;
        res(idx_acq) = 0; % apply another Im
        res = res(idx_nacq);
	end
