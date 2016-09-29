function x = softThresh(y,t)
% apply joint sparsity soft-thresholding

if iscell(y) == 0
    error('y must be a cell array');
end

x = y;

% number of coil
numOfCoil = size(y,1);

% number of wavelet coefficients
numCoeff = numel(y{1,1}.dec);
        
% first coefficient is the low frequency
for k=2:numCoeff
    s = size(y{1,1}.dec{k});
    wCoeff = zeros([s numOfCoil]);
    for n=1:numOfCoil
        wCoeff(:,:,n) = y{n,1}.dec{k} + i*y{n,2}.dec{k};
    end
    
    absy = sqrt(sum(abs(wCoeff).^2,3));
    unity = wCoeff./(repmat(absy,[1,1,numOfCoil])+eps);
    %nity = wCoeff;

    res = absy-t;
    res = (res + abs(res))/2;
    wCoeff = unity.*repmat(res,[1,1,numOfCoil]);
    
    for n=1:numOfCoil
        x{n,1}.dec{k} = real(wCoeff(:,:,n));
        x{n,2}.dec{k} = imag(wCoeff(:,:,n));
    end
end
