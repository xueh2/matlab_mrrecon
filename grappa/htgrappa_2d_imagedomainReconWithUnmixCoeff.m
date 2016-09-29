
function unWarppedCombinedIm = htgrappa_2d_imagedomainReconWithUnmixCoeff(underSampledKSpace, header, unmixCoeff, zeroFilledSize)
%
% unWarppedCombinedIm = htgrappa_2d_imagedomainReconWithUnmixCoeff(underSampledKSpace, header, unmixCoeff)
% This function performs the htgrappa recon using multiplication to replace convolution
% 
% ------------------------ input ------------------------
% underSampledKSpace : undersampled k-space, [Nfe*Npe*numOfCoil], the missing k-space lines are filled with 0 
% header: information header
% header.blocks is the block size for grappaCoeff
% header.Nfe : number of frequency encoding points
% header.Npe : number of full phase encoding lines
% unmixCoeff : image domain unmixing coeffecients
% zeroFilledSize : [newNfe newNpe] if nonempty, zero-filling interpolate the output image
% 
% ------------------------ output ------------------------
% unWarppedCombinedIm : unwarpped complex image with adaptive combination using sensitivity map
% 
% Hui Xue
% March 18, 2011
%
% References:   HTGRAPPA: real-time B1-weighted image domain TGRAPPA reconstruction MRM 61(6) 2009
% ========================================================

numOfCoils = header.num_coils;
Nfe = header.Nfe;
Npe = header.Npe;

s = size(underSampledKSpace);

% the aliased images only need to be computed once
aliasedIm = zeros(s);
for c=1:numOfCoils
    % first, shift 0-N-1
    % do ifft
    % shift back to -N/2 to N/2
    aliasedIm(:,:,c) = ifftshift( ifft2( ifftshift(underSampledKSpace(:,:,c)) ) );
end
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

for cprime = 1:numOfCoils
    unWarppedCombinedIm = unWarppedCombinedIm + aliasedIm(:,:,cprime).*unmixCoeff(:,:,cprime);
end
% figure;imshow(abs(unWarppedCombinedIm(:,:,1)), []);
