function x = nestaSPIRiT(y, GOP, nIter, x0, params)
%-----------------------------------------------------------------------
%
% x = nestaSPIRiT(y, GOP, nIter, lambda, x0)
%
% solve the SPIRiT recon using NESTA algorithm
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
% Phi(x) = ||(G-I)*Dc'*yu+(G-I)*D'*y||^2 + lambda*|W*F'*(Dc'*yu+D'*y)|_1
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

x = x0;

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

mu = 0.01; %--- can be chosen to be small
opts = [];
opts.maxintiter = 5;
opts.TOlVar = 1e-5;
opts.verbose = 0;
opts.maxiter = nIter;
opts.U = @U;
opts.Ut = @Ut;
opts.stoptest = 1;  
opts.typemin = 'l1';
delta = 0;

% testing
coeff = U(y(:));
yu = Ut(coeff);
yu2 = reshape(yu, size(y));
norm(yu2(:)-y(:))

ya = Ac(y(:));
ya2 = Atc(ya);
ya2 = reshape(ya2, size(y));
norm(ya2(:)-y(:))

tic;
[x_nesta,niter,resid,err] = NESTA(@Ac,@Atc,b,mu,delta,opts);
t.NESTA = toc
x = reshape(x_nesta, size(y));
x(idx_acq) = y(idx_acq);

    % -------------------------------------------------------------------
    % define the U and Ut
    function coeff = U(x)
        tempY = reshape(x, size(params.y));
        Im = ifft2c(tempY);
        
        [sx,sy,nc] = size(Im);
        ssx = 2^ceil(log2(sx)); 
        ssy = 2^ceil(log2(sy));
        ss = max(ssx, ssy);
        Im = zpad(Im,ss,ss,nc);
        sparseCoeff = params.sparseTransform*Im;
        coeff = sparseCoeff(:);
    end

    function x = Ut(coeff)
        [sx,sy,nc] = size(params.y);
        ssx = 2^ceil(log2(sx)); 
        ssy = 2^ceil(log2(sy));
        ss = max(ssx, ssy);        
        sparseCoeff = reshape(coeff, [ss ss nc]);
        Im = params.sparseTransform'*sparseCoeff;
        Im = fft2c(crop(Im,sx,sy,nc));
        x = Im(:);
    end

    % define Ac and Atc
    function y2 = Ac(x)
        tempY = reshape(x, size(params.y));
        y2 = reshape(GOP*tempY, [Nfe*Npe*nCoils 1]);
    end

    function x = Atc(y2)
        tempY = reshape(y2, size(params.y));
        x = reshape(GOP'*tempY, [Nfe*Npe*nCoils 1]);
    end

    % -------------------------------------------------------------------

end
