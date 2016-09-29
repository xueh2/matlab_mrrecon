function [m, totalNorm] = coeffNorm(y)
% compute the norm of wavlet coefficients
% using joint sparsity

if iscell(y) == 0
    error('coeff must be a cell array');
end

% number of coil
numOfCoil = size(y,1);

% number of wavelet coefficients
numCoeff = numel(y{1,1}.dec);
s = size(y{1,1}.dec{1});
m = cell([numCoeff-1 1]);

% do not include approximation coefficients
for k=2:numCoeff
    
    s = size(y{1,1}.dec{k});
    m{k-1} = zeros(s);
    tempM = zeros(s);
    
    for n=1:numOfCoil            
        tempM = y{n,1}.dec{k} + i*y{n,2}.dec{k};
        m{k-1} = m{k-1} + tempM.*conj(tempM);
    end
end

totalNorm = 0;
for k=2:numCoeff
    t = sqrt(m{k-1}+1e-15);
    totalNorm = totalNorm + sum(t(:));
end

