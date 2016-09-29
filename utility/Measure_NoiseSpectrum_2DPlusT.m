
function [noiseSpectrum, noiseCorrMatrix] = Measure_NoiseSpectrum_2DPlusT(data, BW, noise, subtractMeanFlag, useNoise)
% This is a function to estimate the noise spectrum from the 2D+T datasets.
% [spectrumNoise, correlationNoise] = Measure_NoiseSpectrum_2DPlusT(data, BW, noise, subtractMeanFlag, useNoise)
% data: data matrix [Nfe Npe frame]
% BW: a mask image to define the region for noise measurement; if empty, the whole image is used
% noise: noise images [Nfe Npe frame], which can be generate using MP law and KLT; if empty, the data matrix is used to generate noise
% subtractMeanFlag: if 1, the mean is removed from data for noise generation

s = size(data);

ind = [];
if ( isempty(BW) )
    BW = data(:,:,1);
    BW(:) = 1;
    ind = find(BW>0);
else
    ind = find(BW>0);
end

N = numel(ind);

COL = s(1);
LIN = s(2);
REP = s(3);

% get the noise
if ( ~isempty(noise) & useNoise )    
    % the noise array is used form compuation
    noiseArray = zeros(N, REP);
    for f=1:REP        
        t = noise(:,:,f);
        noiseArray(:,f) = t(ind(:));        
    end        
else
    noiseArray = zeros(N, REP);
    for f=1:REP        
        t = data(:,:,f);
        noiseArray(:,f) = t(ind(:));        
    end
    
    if ( subtractMeanFlag )
        meanArray = mean(noiseArray, 1);               
        meanArray = repmat(meanArray, [N 1]);                               
        noiseArray = noiseArray - meanArray;
    end        
end

% compute noise spectrum and correlation
noiseSpec = fft(noiseArray, [], 2);        
noiseSpecShifted = fftshift(noiseSpec, 2);
noiseSpecShifted = abs(noiseSpecShifted);      
noiseSpecMean = mean(noiseSpecShifted, 1);
noiseSpectrum = noiseSpecMean;

% noiseCorrMatrix=(1/(M-1))*(noiseArray'*noiseArray);
% noiseCorrMatrix=(noiseCorrMatrix+noiseCorrMatrix')/2;
noiseCorrMatrix = cov(noiseArray);
