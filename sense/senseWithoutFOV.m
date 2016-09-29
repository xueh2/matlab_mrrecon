function [ImCombined, Im, eigD, senMap, x0] = senseWithoutFOV(kspace, ACS, nIter, lambda, senMap, x0)
% perform the iterative sense recon without FOV limitation
% kspace: [COL LIN CHA]
% ACS: fully sampled kspace region for sensitivity estimation

%% estimate the coil sensitivity
imageSize = size(kspace);

if ( nargin < 5 )
    opts=[];
    opts.choice=1; % use opts.tau
    opts.percentage=96;
    kSize=[6,6];

    opts.thd=0.98;
    tic;
    [senMap, eigD, tau]=ED_eigen_2D_parallel_fov(ACS, kSize, imageSize, opts);
    toc;
else
    eigD = [];
end

if ( nargin < 6 )
    % initials
    Im1 = SensitivityCoilCombination(ifft2c(kspace), senMap(:,:,:,1));
    Im2 = SensitivityCoilCombination(ifft2c(kspace), senMap(:,:,:,2));    
    x0 = cat(3, Im1, Im2); 
end

%% perform the recon
% [Im, RESVEC] = cgSENSEWithoutFOV(kspace, senMap, nIter, lambda, x0);

[Im, iter, RESVEC] = cgSENSEWithoutFOV_NotLSQR(kspace, senMap, nIter, lambda, x0);

%% perform final coil combination
% nCoils = size(kspace, 3);
% Im1 = senMap(:,:,:,1) .* repmat(Im(:,:,1), [1 1 nCoils]);
% Im2 = senMap(:,:,:,2) .* repmat(Im(:,:,2), [1 1 nCoils]); 
% ImCombined = Im1 + Im2;
% ImCombined = SensitivityCoilCombination(ImCombined, senMap(:,:,:,1));
ImCombined = Im(:,:,1);

% ImCombined = sum(Im, 3); % plotComplexImage(ImCombined, [1 1 1], 200, 400);
