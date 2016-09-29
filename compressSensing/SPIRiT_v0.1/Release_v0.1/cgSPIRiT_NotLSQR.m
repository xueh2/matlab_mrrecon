function [res, iter] = cgSPIRiT_NotLSQR(y, GOP, nIter, lambda, x0)
%
%
%  res = cgSPIRiT_NotLSQR(y, GOP, nIter, lambda, x0)
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
% The pure conjugate gradient optimization is used here.


tCGSPIRiT=tic;

print_residual = 1;

[sx,sy,nCoils] = size(y);
idx_acq = find(abs(y)>0);
idx_nacq = find(abs(y)==0);

if nargin < 5
     x0 = y;
end

% compute b vector
b = GOP*y; 
b = GOP'*b;
b = -b(idx_nacq(:));
b0 = b;

%% cg iteration
x = x0(idx_nacq(:));
iter = 0;
r = GOP'*(GOP*x0);
r = b - r(idx_nacq(:));
d = r;
delta_new = abs(r'*r);
delta_0 = delta_new;

% bregmanStep = 1;

x_prev = x;
delta_old = delta_new;
v = 0;
while (iter<nIter) & (delta_new>lambda*delta_0) & (delta_new<=delta_old)
    q = y;
    q(:) = 0;
    q(idx_nacq(:)) = d;
    q = GOP'*(GOP*q);
    q = q(idx_nacq(:));
    
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
    
    if (print_residual == 1),
        fprintf('Iteration %d, norm(r) /norm(r0) = %12.8e\n', iter, delta_new/delta_0);drawnow;
    end            
    
    iter = iter + 1;
    
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

res = y;
res(idx_nacq) = x;

disp(['cgSPIRiT_NotLSQR : ' num2str(toc(tCGSPIRiT))]);
