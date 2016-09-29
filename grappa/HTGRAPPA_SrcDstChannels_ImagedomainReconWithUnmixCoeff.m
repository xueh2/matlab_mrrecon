
function unWarppedCombinedIm = HTGRAPPA_SrcDstChannels_ImagedomainReconWithUnmixCoeff(underSampledKSpace, unmixCoeff, zeroFilledSize)
%
% unWarppedCombinedIm = HTGRAPPA_SrcDstChannels_ImagedomainReconWithUnmixCoeff(underSampledKSpace, unmixCoeff, zeroFilledSize)
% This function performs the htgrappa recon using multiplication to replace convolution
% 
% ------------------------ input ------------------------
% underSampledKSpace : undersampled k-space, [Nfe*Npe*srcCha], the missing k-space lines are filled with 0 
% unmixCoeff : image domain unmixing coeffecients, [Nfe Npe srcCha]
% zeroFilledSize : [newNfe newNpe] if nonempty, zero-filling interpolate the output image
% 
% ------------------------ output ------------------------
% unWarppedCombinedIm : unwarpped complex image with adaptive combination using sensitivity map
% 
% Hui Xue
% Oct 10, 2011
%
% References:   HTGRAPPA: real-time B1-weighted image domain TGRAPPA reconstruction MRM 61(6) 2009
% ========================================================

Nfe = size(underSampledKSpace, 1);
Npe = size(underSampledKSpace, 2);
srcCha = size(underSampledKSpace, 3);

if ( nargin < 3 )
    zeroFilledSize = [];
end

s = size(underSampledKSpace);

% the aliased images only need to be computed once
% R = numel(underSampledKSpace) / numel(find(abs(underSampledKSpace)>0));
aliasedIm = ifft2c(underSampledKSpace);
% combinedAliasedIm = SoS_Image(aliasedIm); figure;imshow(abs(combinedAliasedIm), []);

if ( ~isempty(zeroFilledSize) )
    if ( zeroFilledSize(1) ~= s(1) | zeroFilledSize(2)~=s(2) )
        aliasedIm = Zero_Padding_Resize_NoFiltering(aliasedIm, zeroFilledSize(1), zeroFilledSize(2)); 
    end
end

% compute the unaliased images
if ( ~isempty(zeroFilledSize) )
    unWarppedCombinedIm = zeros(zeroFilledSize(1), zeroFilledSize(2));
else
    unWarppedCombinedIm = zeros(Nfe, Npe);
end

% for c = 1:srcCha
%     unWarppedCombinedIm = unWarppedCombinedIm + aliasedIm(:,:,c).*unmixCoeff(:,:,c);
% end
unWarppedCombinedIm = sum(aliasedIm.*unmixCoeff, 3);
% figure;imshow(abs(unWarppedCombinedIm(:,:,1)), []);
