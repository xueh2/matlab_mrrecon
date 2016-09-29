function kernelIm = sgrappaKernelToImage(kernel, imSize, show, voxelsize, centre, width)
% convert the sgrappa kernel into image space using ifft2
% kernel : a kSize(1)*kSize(2)*nCoil*nCoil matrix, corresponding to nCoil sgrappa kernels
% imSize : the image size for recon [row col]
% kernelIm : imSize*nCoil*nCoil, nCoil image space sgrappa kernel
% if show is 1, plot nCoil sgrappa kernel

if nargin < 4
    voxelsize = [ 1 1 1 ];
    centre = 1024;
    width = 1024;
end

kernelIm = zeros([imSize(1) imSize(2) size(kernel,3) size(kernel,4)]);

% numofRow = size(kernel, 1);
% numofCol = size(kernel, 2);
% 
% kernalLengthFE = [-(numofCol-1)/2 : (numofCol+1)/2];
% blockCoord = [-(numofRow-1)/2 : -1];
% blockCoord = [blockCoord 1:(numofRow-1)/2];
% missingLines = 0;
% 
% header.num_coils = size(kernel, 4);
% header.subsampling_factor = 1;
% blocks = header.blocks;
% halfBlock = floor(kernalLengthFE/2);
% Nfe = header.Nfe;
% Npe = header.Npe;

for n=1:size(kernel,4)
    kernelIm(:,:,:,n) = ifft2c(zpad(kernel(end:-1:1,end:-1:1,:,n)*sqrt(imSize(1)*imSize(2)), imSize(1), imSize(2), size(kernel,3)));
end

% kernelIm = htgrappa_2d_imagedomaincoefficients(sgrappaCoeff, header, kernalLengthFE, blockCoord, missingLines);

if show
    for n=1:min(size(kernel,4), 5)
        plotComplexImageArrayMontage(kernelIm(:,:,:,n), voxelsize, centre, width, 1, 1, 1);
        title(['SGRAPPA kernel ' num2str(n)]);
    end
end
