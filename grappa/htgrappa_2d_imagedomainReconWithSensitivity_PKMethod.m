
function unWarppedCombinedIm = htgrappa_2d_imagedomainReconWithSensitivity_PKMethod(underSampledKSpace, sensitivityMap, header, grappaCoeffinIm)
%
% unwarppedIm = htgrappa_2d_imagedomainReconWithSensitivity_PKMethod(underSampledKSpace, sensitivityMap, header, grappaCoeffinIm)
% This function performs the htgrappa recon using multiplication to replace convolution
% The PK method of HTGRAPPA is to perform coil combination on the image domain grappa kernel, then perform the unmixing.
% 
% ------------------------ input ------------------------
% underSampledKSpace : undersampled k-space, [Nfe*Npe*numOfCoil], the missing k-space lines are filled with 0 
% sensitivityMap : [Nfe*Npe*numOfCoil], sensitivity map for each coil 
% header: information header
% header.blocks is the block size for grappaCoeff
% header.Nfe : number of frequency encoding points
% header.Npe : number of full phase encoding lines
% grappaCoeffinIm : grappa kernel in image domain [Nfe*Npe*header.num_coils*header.num_coils]
% e.g. grappaCoeffinIm(:, :, 3, 2) is the 3rd image domain grappa kernel for coil 2
% 
% ------------------------ output ------------------------
% unWarppedCombinedIm : unwarpped complex image with adaptive combination using sensitivity map
% 
% Hui Xue
% March 18, 2011
%
% References:   HTGRAPPA: real-time B1-weighted image domain TGRAPPA reconstruction MRM 61(6) 2009
%               IcePAT slide, Peter Kelleman
% ========================================================

numOfCoils = header.num_coils;
Nfe = header.Nfe;
Npe = header.Npe;

s = size(underSampledKSpace);
scoeff = size(grappaCoeffinIm);

% the aliased images only need to be computed once
aliasedIm = zeros(s);
for c=1:numOfCoils
    % first, shift 0-N-1
    % do ifft
    % shift back to -N/2 to N/2
    aliasedIm(:,:,c) = ifftshift( ifft2( ifftshift(underSampledKSpace(:,:,c)) ) );
end
% combinedAliasedIm = SoS_Image(aliasedIm); figure;imshow(abs(combinedAliasedIm), []);

% compute the per-coil adaptive combined image domain grappa kernel
combinedGrappaCoeffinIm = zeros(s);
for cprime = 1:numOfCoils        
    for c=1:numOfCoils
        combinedGrappaCoeffinIm(:,:,c) = combinedGrappaCoeffinIm(:,:,c) + grappaCoeffinIm(:,:,c, cprime) .* conj(sensitivityMap(:,:,cprime));
    end
end
% figure;imshow(abs(combinedGrappaCoeffinIm(:,:,1)), []);

% compute the unaliased images
unWarppedCombinedIm = zeros(Nfe, Npe);
for cprime = 1:numOfCoils
    unWarppedCombinedIm = unWarppedCombinedIm + aliasedIm(:,:,cprime).*combinedGrappaCoeffinIm(:,:,cprime);
end
% figure;imshow(abs(unWarppedCombinedIm(:,:,1)), []);
