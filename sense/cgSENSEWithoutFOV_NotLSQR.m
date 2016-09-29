function [res, iter, RESVEC] = cgSENSEWithoutFOV_NotLSQR(y, senMap, nIter, lambda, x0)
%
%
%  res = cgSENSEWithoutFOV_NotLSQR(y, senMap, nIter, lambda, x0)
%  
%  Implementation of the Cartesian, conjugate gradiend SENSE reconstruction without FOV limitation
%
%  Input:
%		y -	Undersampled k-space data. Make sure that empty entries are zero
%			or else they will not be filled.
%		senMap - the sensitivity map, including two maps
%		nIter -	Maximum number of iterations
%		lambda-	the stop threshold, here residual norm estimates is RES and the stopping criteria is (RES(n)-RES(n-1))/RES(1) <= lambda
%       x0    - Initil value
%
% Outputs:
%		res - complex images reconed
%       RESVEC -  a vector of the (RES(n)-RES(n-1))/RES(1) at each iteration
% The pure conjugate gradient optimization is used here.


tCGSENSE=tic;

print_residual = 1;

tCGSENSE=tic;

[sx,sy,nCoils] = size(y);

if nargin < 4
	lambda = 1e-4;
end

idx_acq = find(abs(y)>0);
idx_nacq = find(abs(y)==0);

af = numel(y)/numel(idx_acq);

if ( (nargin<5) || isempty(x0) )
    im1 = SensitivityCoilCombination(sqrt(af)*ifft2c(y), senMap(:,:,:,1));
    im2 = SensitivityCoilCombination(sqrt(af)*ifft2c(y), senMap(:,:,:,2));    
    x0 = cat(3, im1, im2);
end

N = length(idx_nacq(:));

b = y(idx_acq(:));
b = AT(b);

%% cg iteration
x = x0(:);
iter = 0;
r = AT(A(x0));
r = b - r;
d = r;
delta_new = abs(r'*r);
delta_0 = delta_new;

% bregmanStep = 1;

RESVEC = [];

x_prev = x;
delta_old = delta_new;
v = 0;
while (iter<nIter)
    q = d;
    q = AT(A(q));
    
    x_prev = x;
    alpha = abs(delta_new/(d'*q));    
    x = x + alpha.*d;
    
    r = r - alpha.*q;
%     currX = parallel.gpu.GPUArray.zeros(size(y));
%     currX(idx_nacq(:)) = x;
%     currX = GOP'*(GOP*currX);
%     r = b - currX(idx_nacq(:)); 
        
    delta_old = delta_new;
    delta_new = abs(r'*r);
    beta = delta_new/delta_old;    
    d = r + beta.*d;    
    
    if ( delta_new > delta_old )
        x = x_prev;
    end
    
    thres = (delta_old-delta_new) / delta_0;
    
    if (print_residual == 1),
        fprintf('Iteration %d, ( norm(r(n))-norm(r(n-1)) ) / norm(r0) = %12.8e\n', iter, thres);drawnow;
    end            
    
    iter = iter + 1;
    
    RESVEC = [RESVEC; thres]; 
    
    if ( abs(thres) <= lambda )
        fprintf('stop ... ');drawnow;
        break;
    end
    
%     if ( mod(iter, bregmanStep) == 0 )
%         % update the residul b
%         currX = parallel.gpu.GPUArray.zeros(size(y));
%         currX(idx_nacq(:)) = x;
%         currR = GOP'*(GOP*currX);
%         v = (b - currR(idx_nacq(:)));
%         b = b0 + v;
%         r = b - currR(idx_nacq(:));
%         d = r;
%         delta_new = abs(r'*r);
%     end
end

res = reshape(x, [sx, sy 2]);

disp(['cgSENSEWithoutFOV_NotLSQR : ' num2str(toc(tCGSENSE))]);

    %% nested function for A and AT
    
    function res = A(x)
        Im = reshape(x, [sx sy 2]);        
        % apply senMap
        Im1 = senMap(:,:,:,1) .* repmat(Im(:,:,1), [1 1 nCoils]);
        Im2 = senMap(:,:,:,2) .* repmat(Im(:,:,2), [1 1 nCoils]);        

        % go to kspace
        res = sqrt(af)*fft2c(Im1+Im2);
        res = res(idx_acq(:));
    end

    function res = AT(x)	
        res = zeros(sx, sy, nCoils);
        res(idx_acq(:)) = x;
        res = sqrt(af)*ifft2c(res);

        Im1 = SensitivityCoilCombination(res, senMap(:,:,:,1));
        Im2 = SensitivityCoilCombination(res, senMap(:,:,:,2));    
        res = cat(3, Im1, Im2); 
        res = res(:);
    end

end
