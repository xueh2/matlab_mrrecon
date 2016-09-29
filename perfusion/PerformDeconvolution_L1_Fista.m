function [r, yr, lambda] = PerformDeconvolution_L1_Fista(cin, y, lambda, useAbsoluteLambda)
% perform the L1 regularized deconvolution with wavelet
% lambda: regularization lamda

%% set up common parameters
N = numel(cin);

if numel(y)~=N
    error('cin and y have different length');
end

if nargin < 4
    useAbsoluteLambda = 0;
end

%% set up the matrix A
A = zeros(N, N);
for i=1:N
    for j=i:-1:1
        A(i,j) = cin(i-j+1);
    end
end

WW=redundantHaarOperator(1);
W = @(x) reshape(WW * reshape(x, [N 1]), 2*N,1);
WT= @(x) reshape(WW'* reshape(x, [2*N 1]), N,1);

opts.n=N;
opts.maxIter=25;
opts.init=0;
if (useAbsoluteLambda)
    opts.rFlag=0; % use the absolute input lamda
    opts.lambda_max = lambda;
else
    opts.rFlag=1;
end
opts.maxFlag=0;
opts.quickWarmStart=0;
opts.scale=5;
opts.verbose = 0;

opts.pFlag=1; % 1: use the one for wavelet
% 0: use the one for L1
opts.W=W;
opts.WT=WT;

if (lambda>0)
    [r, funVal, ValueL] = fista_slep_1D_wavelet(A, A', y, lambda, opts);
else
    % use L curve to find out best lamda
    lambda = -3:0.2:-0.5;
    lambda = 10.^lambda;
    
    resNorm = zeros(numel(lambda), 3);
    
    rAll = zeros(numel(lambda), N);
    
    for ii=numel(lambda):-1:1
        ii;
        [r, funVal, ValueL] = fista_slep_1D_wavelet(A, A', y, lambda(ii), opts);
        rAll(ii, :)  = r(:);
        d = A*r-y;
        
        p = numel(lambda) - ii + 1;
        resNorm(p,1) = lambda(p);
        resNorm(p,2) = norm(d);
        
        d = W(r);
        resNorm(p,3) = norm(d, 1);        
    end
    
    % resNorm(ii,3) = resNorm(ii,3).^2;
    % figure; loglog(resNorm(:,2), resNorm(:,3), '+-');
    
    [k_corner,info] = corner(resNorm(:,2),resNorm(:,3));
    % [reg_c,rho_c,lambda] = l_corner(resNorm(:,2),resNorm(:,3),resNorm(:,1));
    
    lambda = resNorm(k_corner,1);
    [r, funVal, ValueL] = fista_slep_1D_wavelet(A, A', y, lambda, opts);
end

yr = A*r;
% figure; hold on; plot(cin); plot(y, '.'); plot(yr, 'r+'); hold off

% sigma_ratio = sm(:,1) ./ sm(1,1);
% ind = find(sigma_ratio<thres_svd);
% if ( ~empty(ind) )
%     r2 = tgsvd (UU,sm,XX,b,1:ind(end));
% end


