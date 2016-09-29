
function unWarppedCombinedIm = HTGRAPPA_SrcDstChannels_ImagedomainReconWithUnmixCoeff_3D(underSampledKSpace, unmixCoeff)
%
% unWarppedCombinedIm = HTGRAPPA_SrcDstChannels_ImagedomainReconWithUnmixCoeff_3D(underSampledKSpace, unmixCoeff)
% This function performs the htgrappa recon using multiplication to replace convolution
% 
% ------------------------ input ------------------------
% underSampledKSpace : undersampled k-space, [Nfe*Npe*srcCha*Npar], the missing k-space lines are filled with 0 
% unmixCoeff : image domain unmixing coeffecients, [Nfe Npe srcCha Npar]
% 
% ------------------------ output ------------------------
% unWarppedCombinedIm : unwarpped complex image with adaptive combination using sensitivity map
% 
% Hui Xue
% Aug 7, 2013
%
% References:   HTGRAPPA: real-time B1-weighted image domain TGRAPPA reconstruction MRM 61(6) 2009
% ========================================================

if ( size(underSampledKSpace, 4) == 1 )
    aliasedIm = ifft2c(underSampledKSpace);
else
    aliasedIm = ifft3c_Permute(underSampledKSpace);
end
% combinedAliasedIm = SoS_Image(aliasedIm); figure;imshow(abs(combinedAliasedIm), []);

unWarppedCombinedIm = squeeze(sum(aliasedIm.*unmixCoeff, 3));
% figure;imshow(abs(unWarppedCombinedIm(:,:,1,1)), []);
