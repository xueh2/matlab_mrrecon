
function [Cutoff, Variance, SNR, gfactor, BW] = Measure_SNR_3D(data, defineROIFlag, noiseMPorEigenValue, BW, noiseFromBW, voxelsize, centre, width)
% This is a function to estimate the noise variance in 3D image using the MP Law.
% [Cutoff, Variance] = MP_Law_NoiseVariance_3D(data);
% data: original data matrix [Nfe Npe frame]
% defineROIFlag : if 1, the roipoly tool is used to select the signal ROI
% noiseMPorEigenValue : if 1, MP law is used to compute noise; if 2, smallest eigenvalue is used for noise variance
% BW: if set, it is the mask to compute signal
% noiseFromBW : if 1, the noise is also computed from this mask, otherwise, computed form the whole image

if( nargin < 4 )
    BW = [];
end

if( nargin < 5 )
    noiseFromBW = 0; % by default, use the global noise
end

if( nargin < 6 )
    voxelsize = [1 1 1];
    centre = 512;
    width = 1024;
end

s = size(data);
[V, D] = KL_Eigenvalue(data);
E = diag(D); %figure(4), hist(E, 24)
[a,V,D, eigImg] = KL_Eigenimage_NoMean(data);

if ( noiseMPorEigenValue == 1)
    [Cutoff, Variance, ks, beta, p_value, H] = KS_Cutoff_2Steps(E, s(1)*s(2));
    p_value
end     

if ( noiseMPorEigenValue == 2 )
    Variance = min(E);
    Cutoff = 1;    
end

gfactor = std(eigImg(:,:,1:Cutoff),0, 3);

if ( ~defineROIFlag )
    SNR = sqrt(max(E)/s(3))/sqrt(Variance);
    % SNR = sqrt(max(E))/sqrt(Variance);
else
    meanD = mean(data, 3);   
    lowR = 0;
    highR = 4096;
    mag = meanD; 
    mag = normalizeImage2Range(mag, lowR, highR);
    mag = normalizeWindowSetting(mag, centre, width);

    % plotMagImage(abs(meanD), voxelsize, centre, width, 1, 0);    
    if ( isempty(BW) )
        h = figure; imshow(mag, []);
        BW = roipoly;
        close(h);
    end    
    signalD = data.*repmat(BW, [1 1 size(data, 3)]);
    ind = find(abs(signalD)>0);
    meanSignal = mean(data(ind(:)));
    
    % compute noise variance only from heart region
    if ( noiseFromBW )
        [V, D] = KL_Eigenvalue(signalD);
        E = diag(D);
        if ( noiseMPorEigenValue == 1)
            [Cutoff, Variance, ks, beta, p_value, H] = KS_Cutoff_2Steps(E, s(1)*s(2));
            p_value
        end     

        if ( noiseMPorEigenValue == 2 )
            Variance = min(E);
            Cutoff = 1;    
        end
        
        % meanSignal = sqrt(max(E));
    end
    SNR = abs(meanSignal)/sqrt(Variance);    
end