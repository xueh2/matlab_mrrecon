function x = fistaSPIRiT(y, GOP, nIter, x0, params)
%-----------------------------------------------------------------------
%
% x = fistaSPIRiT(y, GOP, nIter, lambda, x0)
%
% solve the SPIRiT recon using FISTA algorithm
%
% The function solves the following problem:
%
% given k-space measurments y, and the SPIRiT operator G
% using the null space concept, the acquired points are y
% the unacquired points are yu
% D is the selection operator for y
% Dc is the selection operator for yu
% W is the wavelet operator -- sparsity transform
% finds the kspace x that minimizes:
%
%
% the optimization method used is non linear conjugate gradient with fast&cheap backtracking
% line-search.
% 
% Hui Xue 2011
%-------------------------------------------------------------------------

if nargin<7
    centre = 1024;
    width = 1024;
end

% data parameters
kernel = getKernel(GOP);
kSize = [size(kernel,1),size(kernel,2)];

[Nfe, Npe, nCoils] = size(y);

numEofY = numel(y);

idx_acq = find(abs(y)>0);
idx_nacq = find(abs(y)==0);
N = length(idx_nacq(:));

% store the acquired data
params.GOP = GOP;
params.y = y;
params.idx_acq = idx_acq;
params.idx_nacq = idx_nacq;
params.nsize = [Nfe Npe nCoils];
  
% do the nesta optimization
b = zeros([numEofY 1]);

STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
stoppingCriterion = STOPPING_OBJECTIVE_VALUE;

maxIteration = params.Itnlim;
lambda = params.xfmWeight;
tolerance = params.objToll;
continuationFlag = 1;
lineSearchFlag = 1;

[sx,sy,nc] = size(y);
ssx = 2^ceil(log2(sx)); 
ssy = 2^ceil(log2(sy));
ss = max(ssx, ssy);

params.kspaceSize = [sx sy nc];
params.xfmSize = [ss ss nc];

imX0 = ifft2c(x0);
imX0 = zpad(imX0,ss,ss,nc);
x0 = params.sparseTransform*imX0;
x0 = reshape(x0, [ss*ss*nc 1]);

clear imX0

tic;
[xp, iterationCount] = SolveFISTA(@Ac, @Atc, b, x0, params, ...
    'maxIteration', maxIteration,...
    'stoppingCriterion', stoppingCriterion,...
    'tolerance', tolerance,...
    'lambda', lambda, ...
    'continuationFlag', continuationFlag, ...
    'lineSearchFlag', lineSearchFlag); 
toc

imR = reshape(xp, [ss ss nc]);
imR = params.sparseTransform'*imR;                        
imR = crop(imR,sx,sy,nc);        
x = fft2c(imR);
x(idx_acq) = y(idx_acq);

    % -------------------------------------------------------------------
    % define Ac and Atc
    function y2 = Ac(x)
        q = reshape(x, [ss ss nc]);
        im = params.sparseTransform'*q;                        
        im = crop(im,sx,sy,nc);        
        im = GOP*fft2c(im);        
        y2 = reshape(im, [Nfe*Npe*nCoils 1]);
    end

    function x = Atc(y2)
        tempY = reshape(y2, size(params.y));
        im = ifft2c(GOP'*tempY);
        im = zpad(im,ss,ss,nc);
        x = params.sparseTransform*im;
        x = reshape(x, [ss*ss*nc 1]);
    end

    % -------------------------------------------------------------------

end
