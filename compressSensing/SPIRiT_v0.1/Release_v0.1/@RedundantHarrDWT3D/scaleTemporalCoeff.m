function coeff = scaleTemporalCoeff(a,coeff,temporalScalingFactor)
% scale the temporal wavelet coefficients

Nrep = size(coeff, 3);
if (Nrep > 1)
    coeff(:,:,Nrep/2+1:Nrep,:) = temporalScalingFactor*coeff(:,:,Nrep/2+1:Nrep,:);
end