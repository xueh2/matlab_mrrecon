function [res, RESVEC] = cgSPIRiT_FullMatrix(y,G, kSize, nIter, lambda, x0)
%
%
%  res = cgSPIRiT_FullMatrix(y,G,nIter, lambda,x0)
%  
%  Implementation of the Cartesian, conjugate gradiend SPIRiT reconstruction
%
%  Input:
%		y -	Undersampled k-space data. Make sure that empty entries are zero
%			or else they will not be filled.
%		G -	the SPIRiT kernel represented by the full size sparse matrix
%		nIter -	Maximum number of iterations
%		lambda-	Tykhonov regularization parameter
%       x0    - Initil value
%
% Outputs:
%		res - Full k-space matrix

tCGSPIRiT=tic;

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

yy = G*y(:)-y(:); 
yy = -yy;

[tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,lambda,nIter, [],[],x0(idx_nacq),G,sx,sy,nCoils,idx_nacq);
% [tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,1e-6,nIter, [],[],x0(idx_nacq),GOP,sx,sy,nCoils,idx_nacq, lambda);
FLAG
RELRES
ITER

res = y;
res(idx_nacq) = tmpres;

disp(['cgSPIRiT_FullMatrix : ' num2str(toc(tCGSPIRiT))]);

function [res,tflag] = aprod(x,G,sx,sy,nCoils,idx_nacq,tflag)
	
	if strcmp(tflag,'transp');
		% tmpy = reshape(x(1:sx*sy*nCoils),sx,sy,nCoils);
        res = G'*x-x;
        res = res(idx_nacq);	
    else
		tmpx = zeros(sx,sy,nCoils);
		tmpx(idx_nacq) = x;
		res = G*tmpx(:)-tmpx(:);
	end
