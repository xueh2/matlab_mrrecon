function [res, RESVEC] = cgSPIRiT3D_GPU(y,GOP, nIter, lambda, x0)
%
%
%  res = cgSPIRiT3D_GPU(y,GOP, nIter, lambda)
%  
%  Implementation of the Cartesian, conjugate gradiend SPIRiT reconstruction
%
%  Input:
%		y -	Undersampled k-space data [COL LIN PAR CHA]. Make sure that empty entries are zero
%			or else they will not be filled.
%		GOP -	the SPIRiT kernel operator obtained by calibration [kfe kpe kpar srcCHA dstCHA]
%		nIter -	Maximum number of iterations
%		lambda-	Tykhonov regularization parameter
%       x0    - Initil value
%
% Outputs:
%		res - Full k-space matrix

tCGSPIRiT=tic;

[sx,sy,sz,nCoils] = size(y);
idx_nacq = find(abs(y)==0);
idx_nacq = gpuArray(single(idx_nacq));

y = gpuArray(y);

if nargin < 4
	lambda = 0;
end

if nargin < 5
     x0 = y;
end
x0 = gpuArray(x0);

yy = GOP*y; 
% yy = [-yy(:); idx_nacq(:)*0];
yy = -yy(:);

[tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,lambda,nIter, [],[],x0(idx_nacq),GOP,sx,sy,sz,nCoils,idx_nacq);

FLAG
RELRES = gather(RELRES)
ITER

res = y;
res(idx_nacq) = tmpres;

res = gather(res);
RESVEC = gather(RESVEC);

disp(['cgSPIRiT3D : ' num2str(toc(tCGSPIRiT))]);

function [res,tflag] = aprod(x,GOP,sx,sy,sz,nCoils,idx_nacq,tflag)
	
	% kernel = getKernel(GOP);
	% kSize = [size(kernel,1),size(kernel,2)];

	if strcmp(tflag,'transp');
		% tmpy = reshape(x(1:sx*sy*nCoils),sx,sy,nCoils);
        res = GOP'*reshape(x(1:sx*sy*sz*nCoils),sx,sy,sz,nCoils);
        % res = res(idx_nacq)+ x(sx*sy*nCoils+1:end)*lambda;
        res = res(idx_nacq);
	
    else
		tmpx = parallel.gpu.GPUArray.zeros(sx,sy,sz,nCoils);
		tmpx(idx_nacq) = x;
		res = GOP*tmpx;
		% res = [res(:) ; lambda*x(:)];
        res = res(:);
	end
