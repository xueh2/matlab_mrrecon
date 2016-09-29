function [R, r, R2] = PerfusionT2StarCorrection_L1(cin_I0, cin_I1, TE0, TE1, lamba)
% R: correction ratio of L1 method
% r: R2* estimation
% R2: correction ratio without L1 constraint

if(nargin<5)
    lamba = 0.005;
end

N = numel(cin_I0)
A = eye(N, N);
A = A*(TE1-TE0);

b = zeros(N, 1);
for f=1:N
    if(cin_I0(f)>cin_I1(f))       
        b(f) = log(cin_I0(f)/cin_I1(f));
    else
        A(f, f) = 0;
    end
end

R2 = zeros(N, 1) + 1;
R2Star = zeros(N, 1) + 0;
for f=1:N
    if(cin_I0(f)>cin_I1(f))       
        R2Star(f) = log(cin_I0(f)/cin_I1(f)) / (TE1-TE0);
        R2(f) = exp(TE0*R2Star(f));
    end    
end

RSNR = cin_I0(:);
RSNR = RSNR.*RSNR.*RSNR.*RSNR;
RSNR = repmat(RSNR, [2 1]);
RSNR_inv = 1.0 ./ RSNR;

WW=redundantHaarOperator(1);
W = @(x) RSNR.*reshape(WW * reshape(x, [N 1]), 2*N,1);
WT= @(x) reshape(WW'* reshape(RSNR_inv.*x, [2*N 1]), N,1);

opts.n=N;
opts.maxIter=100;
opts.init=R2Star;
opts.rFlag=1;
opts.maxFlag=0;
opts.quickWarmStart=0;
opts.scale=5;
opts.keepLL=0;
opts.verbose = 0;
opts.x0 = R2Star;
opts.tol = 1e-9;

opts.pFlag=1; % 1: use the one for wavelet
% 0: use the one for L1
opts.W=W;
opts.WT=WT;

% norm(cin_I0(:))
% norm(cin_I1(:))
% norm(A(:))
% norm(b(:))
% norm(RSNR(:))

[r, funVal, ValueL] = fista_slep_1D_wavelet(A, A', b, lamba, opts);

R = zeros(N, 1);
for f=1:N
    R(f) = exp(TE0*r(f));
end
