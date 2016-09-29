function kernelComposite = computeCompositeNullKernel(kernelIm, show, voxelsize, centre, width)
% convert the pruno composited kernel in image space
% kernelIm : a ImSize(1)*ImSize(2)*nCoil*r matrix, corresponding to r pruno kernels in image space
% kernelComposite : imSize*nCoil*nCoil, nCoil composited image space kernel
% if show is 1, plot r pruno kernel

if nargin < 4
    voxelsize = [ 1 1 1 ];
    centre = 1024;
    width = 1024;
end

Nfe = size(kernelIm,1);
Npe = size(kernelIm,2);
nCoil = size(kernelIm,3);
r = size(kernelIm,4);

kernelComposite = zeros([Nfe Npe nCoil nCoil], 'single');

conjKernelIm = conj(kernelIm);

for i=1:nCoil % composite kernel index
    for j=1:nCoil % subkernel index within main kernel, e.g. kernelComposite(:,:,j,i) is the jth subkernel in ith pruno kernel         
        for k=1:r            
            % kernelComposite(:,:,j,i) = kernelComposite(:,:,j,i) + conj(kernelIm(:,:,i,k)).*kernelIm(:,:,j,k);
            kernelComposite(:,:,j,i) = kernelComposite(:,:,j,i) + conjKernelIm(:,:,i,k).*kernelIm(:,:,j,k);
        end        
    end
end

if show
    for n=1:min(nCoil, 5)
        plotComplexImageArrayMontage(kernelComposite(:,:,:,n), voxelsize, centre, width, 1, 1, 1);
        title(['Composited Kernel ' num2str(n)]);
    end
end