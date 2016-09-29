function res = mtimes(a,x)
% This method performs the SPIRiTSrcDst operator
kernel = a.kernel;
srcCha = size(kernel,3);
dstCha = size(kernel,4);
kSize = size(kernel); 
kSize = kSize(1:2);
method = a.method;

switch method

case {'fft'}
	
	res = zeros(size(x));
    if a.adjoint        
        xx = ifft2c(x);
        for n=1:srcCha
            res(:,:,n) = sum(a.CONJKERNEL(:,:,:,n).*xx,3); 
        end
        res = fft2c(res)-x;
    else
        xx = ifft2c(x);        
        for n=1:dstCha
            res(:,:,n) = sum(a.KERNEL(:,:,:,n).*xx,3);
        end        
        res = fft2c(res)-x;
    end
    
case {'fft_AllGPU'}
	
    fctr = size(x,1)*size(x,2);
   
    if a.adjoint        
        G = sqrt(fctr)*fftshift(ifft2(ifftshift(x)));
        G = repmat(G, [1 1 1 srcCha]);
        G = squeeze(sum(a.CONJKERNELGPU.*G, 3));
        G = 1/sqrt(fctr)*fftshift(fft2(ifftshift(G)));
        res = G-x;
    else
        G = sqrt(fctr)*fftshift(ifft2(ifftshift(x)));
        G = repmat(G, [1 1 1 dstCha]);
        G = squeeze(sum(a.KERNELGPU.*G, 3));
        G = 1/sqrt(fctr)*fftshift(fft2(ifftshift(G)));
        res = G-x;
    end

 case {'image'}
    res = zeros(size(x));
    if a.adjoint
        for n=1:srcCha
            tmpk = squeeze(a.CONJKERNEL(:,:,:,n));
            res(:,:,n) = sum(tmpk.*x,3); 
        end
        res = res - x;
    else
        for n=1:dstCha
            tmpk = a.KERNEL(:,:,:,n);
            res(:,:,n) = sum(tmpk.*x,3); 
        end
        res = res-x;
    end    
end
