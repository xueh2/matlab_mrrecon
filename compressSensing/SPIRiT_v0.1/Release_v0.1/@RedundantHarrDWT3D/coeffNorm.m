function [m, totalNorm] = coeffNorm(y)
% compute the norm of wavlet coefficients
% using joint sparsity

% do not include approximation coefficients
% for k=2:numCoeff   
%     tempM = y(:,:,k,:);
%     m(:,:,k-1) = sum(tempM.*conj(tempM), 4);
% end

m = sum(y.*conj(y), 4);
totalNorm = sum(sqrt(m(:)));
