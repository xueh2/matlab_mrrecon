
function unmixCoeff = htgrappa_2d_computeUnmixingCoeff(sensitivityMap, header, grappaCoeffinIm)
%
% unmixCoeff = htgrappa_2d_computeUnmixingCoeff(sensitivityMap, header, grappaCoeffinIm)
% This function compute the unmixing coefficients for htgrappa recon
% 
% ------------------------ input ------------------------
% sensitivityMap : [Nfe*Npe*numOfCoil], sensitivity map for each coil 
% header: information header
% header.blocks is the block size for grappaCoeff
% header.Nfe : number of frequency encoding points
% header.Npe : number of full phase encoding lines
% grappaCoeffinIm : grappa kernel in image domain [Nfe*Npe*header.num_coils*header.num_coils]
% e.g. grappaCoeffinIm(:, :, 3, 2) is the 3rd image domain grappa kernel for coil 2
% 
% ------------------------ output ------------------------
% unmixCoeff : unmixing coefficients for every coil
% 
% Hui Xue
% March 18, 2011
%
% References:   HTGRAPPA: real-time B1-weighted image domain TGRAPPA reconstruction MRM 61(6) 2009
% ========================================================

numOfCoils = header.num_coils;
Nfe = header.Nfe;
Npe = header.Npe;

s = size(sensitivityMap);
scoeff = size(grappaCoeffinIm);

% unmixing coefficients in image domain
unmixCoeff = zeros(s);
for c=1:numOfCoils
    for cprime = 1:numOfCoils        
        unmixCoeff(:,:,c) = unmixCoeff(:,:,c) + grappaCoeffinIm(:,:,c, cprime) .* conj(sensitivityMap(:,:,cprime));
    end
end
% figure;imshow(abs(unmixCoeff(:,:,1)), []);
