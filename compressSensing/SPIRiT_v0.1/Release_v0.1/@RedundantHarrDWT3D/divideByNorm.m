function coeff = divideByNorm(coeff, coeffNorm, p, mu)
% divided every wavelet coefficients by norm of joint sparsity

numOfCoil = size(coeff, 4);
m = 1./sqrt(coeffNorm + mu);
m = repmat(m, [1 1 1 numOfCoil]);

s = size(coeff);
LLL = coeff(1:s(1)/2, 1:s(2)/2, 1:s(3)/2, :);
coeff = coeff .* m;
coeff(1:s(1)/2, 1:s(2)/2, 1:s(3)/2, :) = LLL;