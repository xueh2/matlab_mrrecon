function [res, RESVEC] = cgSPIRiT_SymmetricLinearEquation_GPU(y,GOP, nIter, lambda, x0)
%
%
%  res = cgSPIRiT_SymmetricLinearEquation_GPU(y,GOP, nIter, lambda)
%  
%  Implementation of the Cartesian, conjugate gradiend SPIRiT reconstruction
%
%  Input:
%		y -	Undersampled k-space data. Make sure that empty entries are zero
%			or else they will not be filled.
%		GOP -	the SPIRiT kernel operator obtained by calibration
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
%   kernel = zeros([kSize,coils,coils]);
%   [AtA,] = corrMatrix(kCalib,kSize);
%   for n=1:coils
%       kernel(:,:,:,n) = calibrate(AtA,kSize,coils,n,0.01);
%   end
%   GOP = SPIRiT(kernel, 'fft',[128,128]);
%
%   % undersample by a factor of 2
%   DATA(1:2:end,2:2:end,:) = 0;
%   DATA(2:2:end,1:2:end,:) = 0;
%   
%   %reconstruct:
%   [res] = cgSPIRiT(DATA,GOP, 20, 1e-5, DATA);
%   figure, imshow(cat(2,sos(imgs), 2*sos(ifft2c(DATA)), sos(ifft2c(res))),[]);
%   title('full,  zero-fill,   result')
%   
% (c) Michael Lustig 2007
%

tCGSPIRiT=tic;

[sx,sy,nCoils] = size(y);
idx_nacq = find(abs(y)==0);
idx_nacq = gpuArray(idx_nacq);

y = gpuArray(y);

if nargin < 4
	lambda = 0;
end

if nargin < 5
     x0 = y;
end
x0 = gpuArray(x0);

yy = GOP'*(GOP*y); 
% yy = [-yy(:); idx_nacq(:)*0];
yy = -yy(idx_nacq);

% [tmpres,FLAG,RELRES,ITER,RESVEC] = lsqrGPU(@aprod,yy,lambda,nIter, [],[],x0(idx_nacq),GOP,sx,sy,nCoils,idx_nacq);
[tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,lambda,nIter, [],[],x0(idx_nacq),GOP,sx,sy,nCoils,idx_nacq);
% [tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,1e-6,nIter, [],[],x0(idx_nacq),GOP,sx,sy,nCoils,idx_nacq, lambda);
FLAG
RELRES = gather(RELRES)
ITER

res = y;
res(idx_nacq) = tmpres;

res = gather(res);
RESVEC = gather(RESVEC);

disp(['cgSPIRiT symmetric equation : ' num2str(toc(tCGSPIRiT))]);

function [res,tflag] = aprod(x,GOP,sx,sy,nCoils,idx_nacq,tflag)
	
	% kernel = getKernel(GOP);
	% kSize = [size(kernel,1),size(kernel,2)];

% 	if strcmp(tflag,'transp');
% 		% tmpy = reshape(x(1:sx*sy*nCoils),sx,sy,nCoils);
%         res = GOP'*reshape(x(1:sx*sy*nCoils),sx,sy,nCoils);
%         res = res(idx_nacq)+ x(sx*sy*nCoils+1:end)*lambda;
%         % res = res(idx_nacq);
% 	
%     else
		tmpx = parallel.gpu.GPUArray.zeros(sx,sy,nCoils);
		tmpx(idx_nacq) = x;
		res = GOP'*(GOP*tmpx);
		% res = [res(:) ; lambda*x(:)];
        res = res(idx_nacq);
% 	end
