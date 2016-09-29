function RawFilterCoeffs = calculateRawFilter (nLINmeas, nCOLmeas, params);
% function RawFilterCoeffs = calculateRawFilter (nLINmeas, nCOLmeas, params);
% 
% function to compute the coefficients for raw filtering (windowing)
% 
% inputs:
%     nLINmeas    number of lines
%     nCOLmeas    number of columns
%     params      structure containing
%         params.type      filter type (e.g., Gaussian, Hanning, Hamming, etc)
%         params.parameter filter parameter (e.g., Gaussian truncation)
% output:
%     RawFilterCoeffs  2d matrix of raw filter coefficients
    
switch params.type
    case {'Gaussian','gaussian'}
        RawFilterCoeffs = gaussian_window([nLINmeas nCOLmeas], params.parameter,'periodic');
    case {'Hanning','hanning'}
        RawFilterCoeffs = hanning(nLINmeas)*hanning(nCOLmeas)';        
    case 'mHanning'
        RawFilterCoeffs = mhanning([nLINmeas nCOLmeas], params.parameter);
    case {'Hamming','hamming'}
        RawFilterCoeffs = hamming(nLINmeas)*hamming(nCOLmeas)';
    case {'Kaiser_Bessel','kaiser_bessel'}
        RawFilterCoeffs = kaiser_bessel_window([nLINmeas nCOLmeas], params.parameter);
    case {'sinus','Sinus'}
        RawFilterCoeffs = sinus([nLINmeas nCOLmeas], params.parameter);
    case 'none'
        RawFilterCoeffs = [];
end

% noise scaling
noise_bandwidth=1/mean(RawFilterCoeffs(:).^2); % noise equivalent Bandwidth
RawFilterCoeffs = RawFilterCoeffs * sqrt(noise_bandwidth);


