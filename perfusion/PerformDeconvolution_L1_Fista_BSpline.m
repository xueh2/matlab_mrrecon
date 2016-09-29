function [r, yr] = PerformDeconvolution_L1_Fista_BSpline(cin, y, orderBSpline, numOfInternalControlPoints, deltaT, lambda, useAbsoluteLambda)
% perform the L1 regularized deconvolution with wavelet and BSpline
% lambda: regularization lamda

%% set up common parameters
N = numel(cin);

if numel(y)~=N
    error('cin and y have different length');
end

if nargin < 7
    useAbsoluteLambda = 0;
end

%% set up BSpline
% number of control points, make sure the clamped bspline is used
if(numOfInternalControlPoints>0)
    p = numOfInternalControlPoints + orderBSpline;
    knots = zeros(1, orderBSpline+p);

    stepSize = (N-1) / (numOfInternalControlPoints+1);

    knots(orderBSpline+1:orderBSpline+numOfInternalControlPoints) = stepSize*(1:numOfInternalControlPoints);
else
    
    M = floor( (N-1)/deltaT );
    
    if ( ((N-1)-M*deltaT)<deltaT/2 )
        M = M -1;
    end
    
    if (M<=1) M = 1; end
    
    numOfInternalControlPoints = M;
    
    knots = zeros(1, orderBSpline+M+orderBSpline);
    
    for kk=1:M
        knots(orderBSpline+kk) = kk*deltaT;
    end
end

p = numOfInternalControlPoints + orderBSpline;
knots(orderBSpline+numOfInternalControlPoints+1:end) = N-1;

B = bspline_basismatrix(orderBSpline, knots, 0:N-1);

%% set up the matrix A
A = zeros(N, N);
for i=1:N
    for j=i:-1:1
        A(i,j) = cin(i-j+1);
    end
end

D = A*B;

WW=redundantHaarOperator(1);
W = @(x) reshape(WW * reshape(x, [p 1]), 2*p,1);
WT= @(x) reshape(WW'* reshape(x, [2*p 1]), p,1);

opts.n=p;
opts.maxIter=20;
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
opts.keepLL = 0;

[coeff, funVal, ValueL] = fista_slep_1D_wavelet(D, D', y, lambda, opts);

yr = D*coeff;
r = B*coeff;

return;

% figure; hold on; plot(cin); plot(y, '.'); plot(yr, 'r+'); hold off

% sigma_ratio = sm(:,1) ./ sm(1,1);
% ind = find(sigma_ratio<thres_svd);
% if ( ~empty(ind) )
%     r2 = tgsvd (UU,sm,XX,b,1:ind(end));
% end

lambda = -4:0.1:1;
lambda = 10.^lambda;

dataNorm = zeros(numel(lambda), 1);
L1Norm = zeros(numel(lambda), 1);

rAll = zeros(numel(lambda), p);

for ii=1:numel(lambda)
    ii
    [r, funVal, ValueL] = fista_slep_1D_wavelet(D, D', y, lambda(ii), opts);
    rAll(ii, :)  = r(:);
    d = D*r-y;
    dataNorm(ii) = norm(d);
    
    d = W(r);
    L1Norm(ii) = norm(d, 1);
    
end
