
function unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff_3D(sensitivityMap, grappaCoeffinIm)
%
% unmixCoeff = HTGRAPPA_SrcDstChannels_UnmixingCoeff_3D(sensitivityMap, grappaCoeffinIm)
% This function compute the unmixing coefficients for htgrappa recon
% 
% ------------------------ input ------------------------
% sensitivityMap : [Nfe*Npe*dstCha*Npar], sensitivity map for each destination channel
% grappaCoeffinIm : grappa kernel in image domain [Nfe*Npe*Npar*srcCha*dstCha]
% e.g. grappaCoeffinIm(:, :, :, 3, 2) is the 3rd image domain grappa kernel for dst coil 2
% 
% ------------------------ output ------------------------
% unmixCoeff : unmixing coefficients for every coil
% 
% Hui Xue
% Aug 7, 2013
%
% References:   HTGRAPPA: real-time B1-weighted image domain TGRAPPA reconstruction MRM 61(6) 2009
% ========================================================

% unmixing coefficients in image domain
Nfe = size(sensitivityMap, 1);
Npe = size(sensitivityMap, 2);
Npar = size(sensitivityMap, 4);
srcCha = size(grappaCoeffinIm, 4);
dstCha = size(grappaCoeffinIm, 5);

unmixCoeff = zeros([Nfe Npe srcCha Npar]);
for c=1:srcCha
    for cprime = 1:dstCha        
        tmp = grappaCoeffinIm(:,:,:,c, cprime) .* squeeze(conj(sensitivityMap(:,:,cprime,:)));
        
        unmixCoeff(:,:,c,:) = unmixCoeff(:,:,c,:) + reshape(tmp, [Nfe Npe 1 Npar]);
    end
end
% figure;imshow(abs(unmixCoeff(:,:,1,1)), []);
