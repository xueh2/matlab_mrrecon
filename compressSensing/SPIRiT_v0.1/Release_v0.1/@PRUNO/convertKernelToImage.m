function kernelIm = convertKernelToImage(kernel, imSize, show, voxelsize, centre, width)
% convert the pruno kernel into image space using ifft2
% kernel : a kSize(1)*kSize(2)*nCoil*r matrix, corresponding to r pruno kernels
% imSize : the image size for recon [row col]
% kernelIm : imSize*nCoil*r, r image space pruno kernel
% if show is 1, plot r pruno kernel

if nargin < 4
    voxelsize = [ 1 1 1 ];
    centre = 1024;
    width = 1024;
end

kernelIm = zeros([imSize(1) imSize(2) size(kernel,3) size(kernel,4)]);

for n=1:size(kernel,4)
    kernelIm(:,:,:,n) = ifft2c(zpad(kernel(end:-1:1,end:-1:1,:,n)*sqrt(imSize(1)*imSize(2)), imSize(1), imSize(2), size(kernel,3)));
end

kernelIm = single(kernelIm);

if show
    for n=1:min(size(kernel,4), 5)
        plotComplexImageArrayMontage(kernelIm(:,:,:,n), voxelsize, centre, width, 1, 1, 1);
        title(['Kernel ' num2str(n)]);
    end
end