function coeff = divideByNorm(coeff, coeffNorm, p, mu)
% divided every wavelet coefficients by norm of joint sparsity
% using joint sparsity

if iscell(coeff) == 0
    error('coeff must be a cell array');
end

% number of coil
numOfCoil = size(coeff,1);

% number of wavelet coefficients
numCoeff = numel(coeff{1,1}.dec);

% do not include the approximation
for k=2:numCoeff
    m = (coeffNorm{k-1} + mu).^(p/2-1);
    for n=1:numOfCoil        
        coeff{n,1}.dec{k} = coeff{n,1}.dec{k} .* m;
        coeff{n,2}.dec{k} = coeff{n,2}.dec{k} .* m;
    end
end
