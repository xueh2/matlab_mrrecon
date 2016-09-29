function [res, RESVEC] = cgSPIRiT2DPlusT(y, GOP, nIter, lambda, x0)
%
%
%  res = cgSPIRiT2DPlusT(y, GOP, nIter, lambda, x0)
%  
%  Implementation of the Cartesian, conjugate gradiend SPIRiT
%  reconstruction for 2D + T
%
%  Input:
%		y -	Undersampled k-space data. Make sure that empty entries are zero
%			or else they will not be filled. [COL LIN CHA REP]
%		GOP -	the SPIRiT kernel operator obtained by calibration, an
%		object of SPIRiT2DPlusT
%		nIter -	Maximum number of iterations
%		lambda-	Tykhonov regularization parameter
%       x0    - Initil value
%
% Outputs:
%		res - Full k-space matrix [COL LIN CHA REP]

tCGSPIRiT=tic;

if nargin < 4
	lambda = 0;
end

if nargin < 5
     x0 = y;
end

[sx,sy,nCoils, nReps] = size(y);

idx_nacq = find(abs(y)==0);

yy = GOP*y; 
% yy = [-yy(:); idx_nacq(:)*0];
yy = [-yy(:)];

[tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,lambda,nIter, [],[],x0(idx_nacq),GOP,sx,sy,nCoils,nReps,idx_nacq);
% [tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,1e-6,nIter, [],[],x0(idx_nacq),GOP,sx,sy,nCoils,idx_nacq, lambda);
FLAG
RELRES
ITER

res = y;
res(idx_nacq) = tmpres;

disp(['cgSPIRiT : ' num2str(toc(tCGSPIRiT))]);

function [res,tflag] = aprod(x,GOP,sx,sy,nCoils,nReps,idx_nacq,tflag)
	
	% kernel = getKernel(GOP);
	% kSize = [size(kernel,1),size(kernel,2)];

	if strcmp(tflag,'transp');
		% tmpy = reshape(x(1:sx*sy*nCoils),sx,sy,nCoils);
        res = GOP'*reshape(x(1:sx*sy*nCoils*nReps),sx,sy,nCoils,nReps);
        % res = res(idx_nacq)+ x(sx*sy*nCoils+1:end)*lambda;
        res = res(idx_nacq);
	
    else
		tmpx = zeros(sx,sy,nCoils,nReps);
		tmpx(idx_nacq) = x;
		res = GOP*tmpx;
		% res = [res(:) ; lambda*x(:)];
        res = res(:);
	end
