function [m, totalNorm] = coeffNorm(y)
% compute the norm of wavlet coefficients
% using joint sparsity

% do not include approximation coefficients
% for k=2:numCoeff   
%     tempM = y(:,:,k,:);
%     m(:,:,k-1) = sum(tempM.*conj(tempM), 4);
% end

m = y(:,:,2:end,:);
m = sum(m.*conj(m), 4);

% totalNorm = sqrt(m);
% totalNorm = totalNorm(:)'*totalNorm(:);
% totalNorm

totalNorm = sum(sqrt(m(:)));

% number of coil
% numOfCoil = size(y,4);
% 
% % number of wavelet coefficients
% numCoeff = size(y,3);
% s = [size(y, 1) size(y, 2)];
% m2 = cell([numCoeff-1 1]);
% 
% % do not include approximation coefficients
% for k=2:numCoeff
%     
%     m2{k-1} = zeros(s);
%     tempM = zeros(s);
%     
%     for n=1:numOfCoil            
%         tempM = y(:,:,k,n);
%         m2{k-1} = m2{k-1} + tempM.*conj(tempM);
%     end
% end
% 
% totalNorm = 0;
% for k=2:numCoeff
%     t = sqrt(m2{k-1}+1e-15);
%     totalNorm = totalNorm + sum(t(:));
% end
% totalNorm
