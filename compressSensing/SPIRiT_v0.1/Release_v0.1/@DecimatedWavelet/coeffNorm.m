function m = coeffNorm(y)
% compute the norm of wavlet coefficients

if iscell(y) == 0
    error('coeff must be a cell array');
end

% number of coil
numOfCoil = size(y,1);
m = zeros(numOfCoil,1);

% number of wavelet coefficients
numCoeff = numel(y{1,1}.dec);

s = size(y{1,1}.dec{1});
wCoeff = y{1,1}.dec{1} + i * y{1,2}.dec{1};
for n=1:numOfCoil
    for k=1:numCoeff        
        wCoeff = y{n,1}.dec{k} + i*y{n,2}.dec{k};
        m(n) = m(n) + wCoeff(:)'*wCoeff(:);
    end
end
