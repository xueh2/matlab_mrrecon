function thres = EstimateSureShrinkThreshold(waveletCoeff, noiseSigma)
% estimate the Sure Threshold for the input wavelet coefficients
% following the paper "A Wavelet Tour of Signal Processing"

N = numel(waveletCoeff);
coeff = abs(waveletCoeff);
% compute the risks
risks = zeros(N, 1);

% "A Wavelet Tour of Signal Processing"
coeffSorted = sort(coeff(:), 1, 'descend');
coeff2 = coeffSorted .* coeffSorted;
risks(1) = sum(coeff2(:)) - (N-1)*noiseSigma; % (noiseSigma*noiseSigma+coeff(1)*coeff(1));
for l=2:N
    risks(l) = risks(l-1) - coeff2(l-1) + noiseSigma; 
end

risks = risks + coeff2(:) + noiseSigma*noiseSigma;

% "Dohono paper"
% coeffSorted = sort(coeff(:));
% coeff2 = coeffSorted .* coeffSorted;
% risks(1) = N - 2*N + sum(coeff2(:));
% for l=2:N
%     risks(l) = risks(l-1) + 2 - coeff2(l); 
% end

% find the threshold
[thres, ind] = min(risks);
thres = coeffSorted(ind);

% in case of signal energy is too small
% episilon = noiseSigma*noiseSigma*sqrt(N)*((log(N))^1.5);
% if ( sum(coeff2(:))-N*noiseSigma <= episilon )
%     thres = noiseSigma*sqrt(2*log(N));
% end

