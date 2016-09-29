
function [combinedIm, sensitivityMap, eigD] = Spirit2DForSensitivity(complexImage, kSize, thres)
% ------------------------ input ------------------------
% complexImage : complex images, [Nfe Npe num_coils)]
%
% ------------------------ output ------------------------
% combinedIm : combined complex image
% sensitivityMap : estimated sensitivity map for each point
% 
% Hui Xue
% Dec 15, 2011
%
% ========================================================
s = size(complexImage);
numOfCoils = s(3);
Nfe = s(1);
Npe = s(2);

sensitivityMap = zeros([s 2]);

kernel = spiritCalibration(fft2c(complexImage), kSize, thres);
GOP = SPIRiT(kernel, 'fft',[Nfe,Npe]);

eigD = zeros(Nfe, Npe, 2);
kIm = GOP.KERNEL;

for pe = 1:Npe
    for fe = 1:Nfe               
        R = reshape(kIm(fe, pe, :, :), [numOfCoils numOfCoils]);        
        [V, D] = eig((R+R')/2);        
        sensitivityMap(fe, pe, :, 1) = V(:, end);
        sensitivityMap(fe, pe, :, 2) = V(:, end-1);
        D = diag(D);
        eigD(fe, pe, 1) = D(end);
        eigD(fe, pe, 2) = D(end-1);
    end
end
sensitivityMap = conj(sensitivityMap);

combinedIm = SensitivityCoilCombination(complexImage, sensitivityMap);
