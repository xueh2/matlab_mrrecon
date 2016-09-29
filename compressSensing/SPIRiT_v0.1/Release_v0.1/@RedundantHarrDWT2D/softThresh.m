function x = softThresh(a,y,t,gFactor)
% apply joint sparsity soft-thresholding

x = y;

% number of coil
numOfCoil = size(y,4);

% number of wavelet coefficients
numCoeff = size(y,3);
        
if ( isempty(gFactor) )
    % first coefficient is the low frequency
    for k=2:numCoeff

        wCoeff = squeeze(y(:,:,k,:));

        absy = sqrt(sum(abs(wCoeff).^2,3));
        unity = wCoeff./(repmat(absy,[1,1,numOfCoil])+eps);

        res = absy-t;
        res = (res + abs(res))/2;
        wCoeff = unity.*repmat(res,[1,1,numOfCoil]);

        x(:,:,k,:) = wCoeff;
    end
else
    % first coefficient is the low frequency
    for k=2:numCoeff

        wCoeff = squeeze(y(:,:,k,:));

        absy = sqrt(sum(abs(wCoeff).^2,3));
        unity = wCoeff./(repmat(absy,[1,1,numOfCoil])+eps);

        tUsed = zeros(size(absy))+t;
        tUsed = tUsed .* gFactor;
        res = absy-tUsed;
        res = (res + abs(res))/2;
        wCoeff = unity.*repmat(res,[1,1,numOfCoil]);

        x(:,:,k,:) = wCoeff;
    end
end