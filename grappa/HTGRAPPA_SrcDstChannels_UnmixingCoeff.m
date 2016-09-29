
function unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff(sensitivityMap, grappaCoeffinIm)
%
% unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff(sensitivityMap, grappaCoeffinIm)
% This function compute the unmixing coefficients for htgrappa recon
% 
% ------------------------ input ------------------------
% sensitivityMap : [Nfe*Npe*dstCha], sensitivity map for each destination channel
% grappaCoeffinIm : grappa kernel in image domain [Nfe*Npe*srcCha*dstCha]
% e.g. grappaCoeffinIm(:, :, 3, 2) is the 3rd image domain grappa kernel for dst coil 2
% 
% ------------------------ output ------------------------
% unmixCoeff : unmixing coefficients for every coil
% 
% Hui Xue
% Oct 10, 2011
%
% References:   HTGRAPPA: real-time B1-weighted image domain TGRAPPA reconstruction MRM 61(6) 2009
% ========================================================

% unmixing coefficients in image domain
Nfe = size(sensitivityMap, 1);
Npe = size(sensitivityMap, 2);
srcCha = size(grappaCoeffinIm, 3);
dstCha = size(grappaCoeffinIm, 4);

unmixCoeff = zeros([Nfe Npe srcCha]);
for c=1:srcCha
    for cprime = 1:dstCha        
        unmixCoeff(:,:,c) = unmixCoeff(:,:,c) + grappaCoeffinIm(:,:,c, cprime) .* conj(sensitivityMap(:,:,cprime));
    end
end
% figure;imshow(abs(unmixCoeff(:,:,1)), []);
