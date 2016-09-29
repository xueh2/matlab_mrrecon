function [res, RESVEC] = cgSPIRiT3D(y,GOP, nIter, lambda, x0)
%
%
%  res = cgSPIRiT3D(y,GOP, nIter, lambda)
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

tCGSPIRiT=tic;

if nargin < 4
	lambda = 0;
end

if nargin < 5
     x0 = y;
end

[sx,sy,sz,nCoils] = size(y);

idx_acq = find(abs(y)>0);
idx_nacq = find(abs(y)==0);
N = length(idx_nacq(:));

yy = GOP*y; 
yy = [-yy(:); idx_nacq(:)*0];
% yy = [-yy(:)];

[tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,1e-6,nIter, [],[],x0(idx_nacq),GOP,sx,sy,sz,nCoils,idx_nacq, lambda);

FLAG
RELRES
ITER

res = y;
res(idx_nacq) = tmpres;

disp(['cgSPIRiT3D : ' num2str(toc(tCGSPIRiT))]);

function [res,tflag] = aprod(x,GOP,sx,sy,sz,nCoils,idx_nacq,lambda,tflag)
	
	% kernel = getKernel(GOP);
	% kSize = [size(kernel,1),size(kernel,2)];

	if strcmp(tflag,'transp');
		% tmpy = reshape(x(1:sx*sy*nCoils),sx,sy,nCoils);
        
        res = GOP'*reshape(x(1:sx*sy*sz*nCoils),sx,sy,sz,nCoils);
        res = res(idx_nacq)+ x(sx*sy*sz*nCoils+1:end)*lambda;
        % res = res(idx_nacq);
	
    else
		tmpx = zeros(sx,sy,sz,nCoils,'single');
		tmpx(idx_nacq) = x;
		res = GOP*tmpx;
		res = [res(:) ; lambda*x(:)];
        % res = res(:);
	end
