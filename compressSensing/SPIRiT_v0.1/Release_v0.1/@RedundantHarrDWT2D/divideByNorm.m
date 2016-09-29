function coeff = divideByNorm(coeff, coeffNorm, p, mu)
% divided every wavelet coefficients by norm of joint sparsity
% using joint sparsity
numOfCoil = size(coeff, 4);
m = 1./sqrt(coeffNorm + mu);
m = repmat(m, [1 1 1 numOfCoil]);
coeff(:,:,2:end,:) = coeff(:,:,2:end,:) .* m;