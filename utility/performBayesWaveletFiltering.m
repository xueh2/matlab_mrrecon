function [ori_KLT, ori_KLT_WavEigen, ori_KLT_OnlyWavEigen, thres, eigImg] = performBayesWaveletFiltering(data, header, numOfWavLevels, wname, ratioForThres, plotFlag)
% N : number of wavelet transform levels
% wname : wavelet name, following the wavelet toolbox
% ratioForThres : ratio to scale the estiamted thresholds, realThresUsed = thres*ratioForThres
% ori_KLT : only temporal KLT
% ori_KLT_WavEigen : temporal KLT + spatial wavelet on eigenImage
% ori_KLT_OnlyWavEigen : only spatial wavelet on eigenImage

[a,V,D, eigImg] = KL_Eigenimage_NoMean(data);

% perform eigen image filtering, first eigen image does not filer
filteredEigIm = a;
thresAll = zeros(header.sizeZ, 3*numOfWavLevels+1);
for eigIm=header.sizeZ:-1:1
    im = a(:,:,eigIm); % imtool(im, [])
    [thres, filteredWaveCoeff, filteredImg] = estimateBayesWaveletThreshold(im, numOfWavLevels, wname, ratioForThres); % imtool(filteredImg, [])
    filteredEigIm(:,:,eigIm) = filteredImg;
    if ( ~isempty(thres) ) 
        thresAll(eigIm, :) = thres(:);
    end
end

eigV = sqrt(diag(D));

if ( plotFlag )
    figure;
    for pp=1:size(thresAll,2)
        subplot(2,ceil(size(thresAll,2)/2),pp);
        hold on
        plot(thresAll(:,pp)./eigV, 'r');
        plot(thresAll(:,pp)./eigV, 'b+');
        hold off
        axis on
        grid on
        box on
    end
end

[ Cutoff2, Variance1, ks, beta_ratio, p_value ] = KS_Cutoff_Eigenimages_2Steps(data);
Cutoff2

Cutoff = Cutoff2

if ( Cutoff<0.1*header.sizeZ | Cutoff>0.9*header.sizeZ )
    Cutoff = round(header.sizeZ*0.5);
end

N = header.sizeZ;
r = Cutoff/N;
V2 = V;
realCutOff = round(N*r);
V2(:,1:round(N*r)) = 0;

s = size(a);
b = (reshape(a, [s(1)*s(2),s(3)] ));
ori_KLT = (reshape(b*V2', [s(1),s(2),s(3)] ) ) ;
% ori_KLT(find(ori_KLT<0)) = 0;

b = (reshape(filteredEigIm, [s(1)*s(2),s(3)] ));
ori_KLT_WavEigen = (reshape(b*V2', [s(1),s(2),s(3)] ) ) ;
% ori_KLT_WavEigen(find(ori_KLT_WavEigen<0)) = 0;

b = (reshape(filteredEigIm, [s(1)*s(2),s(3)] ));
ori_KLT_OnlyWavEigen = (reshape(b*V', [s(1),s(2),s(3)] ) ) ;
% ori_KLT_OnlyWavEigen(find(ori_KLT_OnlyWavEigen<0)) = 0;
