
function unWarppedIm = htgrappa_2d_imagedomainRecon_InputImage(aliasedIm, header, grappaCoeffinIm)
%
% unwarppedIm = htgrappa_2d_imagedomainRecon_InputImage(aliasedIm, header, grappaCoeffinIm)
% This function performs the htgrappa recon using multiplication to replace convolution
% 
% ------------------------ input ------------------------
% aliasedIm : aliased image [Nfe*Npe*numOfCoil]
% header: information header
% header.blocks is the block size for grappaCoeff
% header.Nfe : number of frequency encoding points
% header.Npe : number of full phase encoding lines
% grappaCoeffinIm : grappa kernel in image domain [Nfe*Npe*header.num_coils*header.num_coils]
% e.g. grappaCoeffinIm(:, :, 3, 2) is the 3rd image domain grappa kernel for coil 2
% 
% ------------------------ output ------------------------
% unWarppedIm : unwarpped complex image for every coil [Nfe*Npe*numOfCoil]
% 
% Hui Xue
% Jan 24, 2011
%
% References:   HTGRAPPA: real-time B1-weighted image domain TGRAPPA reconstruction MRM 61(6) 2009
% ========================================================

numOfCoils = header.num_coils;
Nfe = header.Nfe;
Npe = header.Npe;

s = size(aliasedIm);
scoeff = size(grappaCoeffinIm);

unWarppedIm = zeros(Nfe, Npe, numOfCoils);

% % the aliased images only need to be computed once
% aliasedIm = zeros(s);
% for c=1:numOfCoils
%     % first, shift 0-N-1
%     % do ifft
%     % shift back to -N/2 to N/2
%     aliasedIm(:,:,c) = ifftshift( ifft2( ifftshift(underSampledKSpace(:,:,c)) ) );
% end
    
% compute the unaliased images
for cprime = 1:numOfCoils
    unaliasedIm = zeros(Nfe, Npe);
    for c=1:numOfCoils
        % suppose the grappaCoeffinIm is -N/2 to N/2
        unaliasedIm = unaliasedIm + aliasedIm(:, :, c) .* grappaCoeffinIm(:,:,c,cprime);
    end
    unWarppedIm(:,:,cprime) = unaliasedIm;
end
