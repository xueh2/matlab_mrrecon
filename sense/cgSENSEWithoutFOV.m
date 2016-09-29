function [res, RESVEC] = cgSENSEWithoutFOV(y, senMap, nIter, lambda, x0)
%
%  res = cgSENSEWithoutFOV(y, senMap, nIter, lambda, x0)
%  
%  Implementation of the Cartesian, conjugate gradiend SENSE reconstruction
%  Without FOV limitation SENSE
%
%  Hui Xue
%  hui.xue@nih.gov
% -------------------------------------------------------------------------

tCGSENSE=tic;

[sx,sy,nCoils] = size(y);

if nargin < 4
	lambda = 0;
end

if nargin < 5    
    im1 = SensitivityCoilCombination(ifft2c(y), senMap(:,:,:,1));
    im2 = SensitivityCoilCombination(ifft2c(y), senMap(:,:,:,2));    
    x0 = cat(3, im1, im2);
end

idx_acq = find(abs(y)>0);
idx_nacq = find(abs(y)==0);
N = length(idx_nacq(:));

b = y(idx_acq(:));

[tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,b,lambda,nIter, [],[],x0(:),senMap,sx,sy,nCoils,idx_acq,idx_nacq);
% [tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,1e-6,nIter, [],[],x0(idx_nacq),GOP,sx,sy,nCoils,idx_nacq, lambda);
FLAG
RELRES
ITER

res = reshape(tmpres, [sx, sy 2]);

disp(['cgSENSEWithoutFOV : ' num2str(toc(tCGSENSE))]);

function [res,tflag] = aprod(x,senMap,sx,sy,nCoils,idx_acq,idx_nacq,tflag)
	
	if strcmp(tflag,'transp');

        res = zeros(sx, sy, nCoils);
        res(idx_acq(:)) = x;
        res = ifft2c(res);
        
        Im1 = SensitivityCoilCombination(res, senMap(:,:,:,1));
        Im2 = SensitivityCoilCombination(res, senMap(:,:,:,2));    
        res = cat(3, Im1, Im2); 
        res = res(:);
    else        
		Im = reshape(x, [sx sy 2]);        
        % apply senMap
        Im1 = senMap(:,:,:,1) .* repmat(Im(:,:,1), [1 1 nCoils]);
        Im2 = senMap(:,:,:,2) .* repmat(Im(:,:,2), [1 1 nCoils]);        
        
        % go to kspace
        res = fft2c(Im1+Im2);
        res = res(idx_acq(:));
	end
