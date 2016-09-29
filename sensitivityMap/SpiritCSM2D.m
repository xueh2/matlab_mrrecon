function [sensitivityMap, eigD] = SpiritCSM2D(kCalib, kSize, imageSize, thres)
% ------------------------ input ------------------------
% kCalib : calibration data,         [cx, cy, num_coils)]
% kSize:   calibration kernel size,  [kx, ky]
% thres:   thre for computing the calibration kernel
%
% ------------------------ output ------------------------
% sensitivityMap : estimated sensitivity map for each point
% 
% Jun Liu
% January 27, 2012
%
% ========================================================

coils=size(kCalib,3);
Nfe=imageSize(1);
Npe=imageSize(2);
% for calibration data: fe, pe

% CalibTyk = 0.01;  % Tykhonov regularization in the calibration
% [AtA,] = corrMatrix(kCalib,kSize);
% for n=1:coils
% 	kernel(:,:,:,n) = calibrate(AtA,kSize,coils,n,CalibTyk);
% end

kernel = spiritCalibration(kCalib, kSize, thres);
GOP = SPIRiT(kernel, 'fft',[Nfe,Npe]);

kIm = GOP.KERNEL;

sensitivityMap = zeros(Nfe, Npe, coils, 2);
eigD = zeros(Nfe, Npe, coils);

for pe = 1:Npe
    for fe = 1:Nfe               
        R = reshape(kIm(fe, pe, :, :), [coils coils]);  
        temp=R-eye(size(R));
        tt = temp * temp';
        [V, D] = eig( tt );  % get the samllest one
        
        sensitivityMap(fe, pe, :, 1) = V(:, 1);
        sensitivityMap(fe, pe, :, 2) = V(:, 2);
        
        D = diag(D);
        
        eigD(fe, pe, :) = D(:);
    end
end
sensitivityMap = conj(sensitivityMap);

