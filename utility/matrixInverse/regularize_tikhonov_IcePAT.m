
function regA = regularize_tikhonov_IcePAT(A, varargin)
% first order tikhonov regularization using IcePAT strategy

s = size(A);

if nargin == 2
    sigma = varargin{1};
else
    sigma = 0.0001; % 1/10000 of the maximum singular value
end

if sigma > 0.99999999,
    'Error! sigma must be a real number (0,1)!'
    return
end
if s(1) ~=s(2)
    'Error!, c_0 must be a square matrix'
    return
end

dTraceR = trace(real(A));
nRow = size(A, 1);
lamda  = sigma * dTraceR / nRow;
for n=1:nRow
    A(n,n) = real(A(n,n)) + lamda;
end
regA = A;
