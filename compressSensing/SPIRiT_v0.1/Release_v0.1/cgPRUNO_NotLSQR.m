function [res, iter] = cgPRUNO_NotLSQR(y, PRU, nIter, lambda, x0)
%
%
%  res = cgPRUNO_NotLSQR(y, GOP, nIter, lambda, x0)
%  
%  Implementation of the Cartesian, conjugate gradiend PRUNO reconstruction
%
%  Input:
%		y -	Undersampled k-space data. Make sure that empty entries are zero
%			or else they will not be filled.
%		PRU -	the PRUNO kernel operator obtained by calibration, performing N'N
%		nIter -	Maximum number of iterations
%		lambda-	Tykhonov regularization parameter
%       x0    - Initil value
%
% Outputs:
%		res - Full k-space matrix
%
% The pure conjugate gradient optimization is used here.


tCGPRUNO=tic;

print_residual = 1;

[sx,sy,nCoils] = size(y);
idx_acq = find(abs(y)>0);
idx_nacq = find(abs(y)==0);

if nargin < 5
     x0 = y;
end

% compute b vector
b = PRU*y; % (N'N)*D'*y 
b(idx_acq) = 0;
b = -b(idx_nacq);

%% cg iteration
x = x0(idx_nacq(:));
iter = 0;
r = PRU*x0;
r = b - r(idx_nacq(:));
d = r;
delta_new = abs(r'*r);
delta_0 = delta_new;

delta_old = delta_new;
while (iter<nIter) & (delta_new>lambda*delta_0) & (delta_new<=delta_old)
    q = y;
    q(:) = 0;
    q(idx_nacq(:)) = d;
    q = PRU*q;
    q = q(idx_nacq(:));
    
    x_prev = x;
    alpha = abs(delta_new/(d'*q));    
    x = x + alpha.*d;
    
    r = r - alpha.*q;
       
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
end

res = y;
res(idx_nacq) = x;
res = single(res);

disp(['cgPRUNO_NotLSQR : ' num2str(toc(tCGPRUNO))]);
