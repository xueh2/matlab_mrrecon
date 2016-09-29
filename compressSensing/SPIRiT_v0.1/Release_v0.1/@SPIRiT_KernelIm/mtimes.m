function res = mtimes(a,x)
% This method performs the SPIRiT operator
kernel = a.kernel;
nCoils = size(x,3);
kSize = size(kernel); kSize = kSize(1:2);
method = a.method;
% KERNEL = a.KERNEL;

switch method

case {'conv'}
	res = zeros(size(x));
    if a.adjoint
        for n=1:nCoils
            tmpk = kernel(:,:,:,n);
            tmpk(floor(kSize(1)/2)+1, floor(kSize(2)/2)+1,n) = ...
			tmpk(floor(kSize(1)/2)+1, floor(kSize(2)/2)+1,n) -1;
            res = res + adjSPIRiT(x(:,:,n),tmpk); 
        end
    else
        for n=1:nCoils
            tmpk = kernel(:,:,:,n);
            tmpk(floor(kSize(1)/2)+1, floor(kSize(2)/2)+1,n) = ...
		    tmpk(floor(kSize(1)/2)+1, floor(kSize(2)/2)+1,n) -1;
            res(:,:,n) = SPIRiT(x, tmpk);
        end

    end

case {'fft'}
	
	res = zeros(size(x));
    if a.adjoint
        
        xx = ifft2c(x);
        for n=1:nCoils
            % res(:,:,n) = sum(squeeze(conj(KERNEL(:,:,n,:))).*xx,3); 
            res(:,:,n) = sum(a.CONJKERNEL(:,:,:,n).*xx,3); 
        end
        res = fft2c(res);
    else
        xx = ifft2c(x);
        
        for n=1:nCoils
            res(:,:,n) = sum(a.KERNEL(:,:,:,n).*xx,3);
        end
        
        res = fft2c(res);
    end
    
case {'fft_MrFtk'}
	
	res = zeros(size(x));
    if a.adjoint        
        xx = ifft2c_MrFtk(x);
        res = Matlab_PerformImageDomainUnwarpping(a.CONJKERNEL, xx);
        res = fft2c_MrFtk(res)-x;
    else
        xx = ifft2c_MrFtk(x);
        res = Matlab_PerformImageDomainUnwarpping(a.KERNEL, xx);
        res = fft2c_MrFtk(res)-x;
    end
    
case {'fft_GPU'}
	
	res = zeros(size(x));
    fctr = size(x,1)*size(x,2);
%     tic; 
    
    % load x into GPU
    G = gpuArray(x);

    if a.adjoint        
%         tic; G = sqrt(fctr)*fftshift(ifft2(ifftshift(G)));toc
%         tic; G = repmat(G, [1 1 1 nCoils]);toc
%         tic; G = squeeze(sum(a.CONJKERNELGPU.*G, 3));toc
%         tic; G = 1/sqrt(fctr)*fftshift(fft2(ifftshift(G)));toc
%         tic; res = gather(G)-x;toc
        G = sqrt(fctr)*fftshift(ifft2(ifftshift(G)));
        G = repmat(G, [1 1 1 nCoils]);
        G = squeeze(sum(a.CONJKERNELGPU.*G, 3));
        G = 1/sqrt(fctr)*fftshift(fft2(ifftshift(G)));
        % res = gather(G)-x;
        res = gather(G);
    else
%         tic; G = sqrt(fctr)*fftshift(ifft2(ifftshift(G)));toc
%         tic; G = repmat(G, [1 1 1 nCoils]);toc
%         tic; G = squeeze(sum(a.KERNELGPU.*G, 3));toc
%         tic; G = 1/sqrt(fctr)*fftshift(fft2(ifftshift(G)));toc
%         tic; res = gather(G)-x;toc
        G = sqrt(fctr)*fftshift(ifft2(ifftshift(G)));
        G = repmat(G, [1 1 1 nCoils]);
        G = squeeze(sum(a.KERNELGPU.*G, 3));
        G = 1/sqrt(fctr)*fftshift(fft2(ifftshift(G)));
        % res = gather(G)-x;
        res = gather(G);
    end

case {'fft_AllGPU'}
	
    fctr = size(x,1)*size(x,2);
   
    if a.adjoint        
        G = sqrt(fctr)*fftshift(ifft2(ifftshift(x)));
        G = repmat(G, [1 1 1 nCoils]);
        G = squeeze(sum(a.CONJKERNELGPU.*G, 3));
        G = 1/sqrt(fctr)*fftshift(fft2(ifftshift(G)));
        res = G;
    else
        G = sqrt(fctr)*fftshift(ifft2(ifftshift(x)));
        G = repmat(G, [1 1 1 nCoils]);
        G = squeeze(sum(a.KERNELGPU.*G, 3));
        G = 1/sqrt(fctr)*fftshift(fft2(ifftshift(G)));
        res = G;
    end
    
 case {'fft_noshift_AllGPU'}
     
     fctr = size(x,1)*size(x,2);
   
    if a.adjoint        
        G = sqrt(fctr)*ifft2(x);
        G = repmat(G, [1 1 1 nCoils]);
        G = squeeze(sum(a.CONJKERNELGPU.*G, 3));
        G = 1/sqrt(fctr)*fft2(G);
        res = G-x;
    else
        G = sqrt(fctr)*ifft2(x);
        G = repmat(G, [1 1 1 nCoils]);
        G = squeeze(sum(a.KERNELGPU.*G, 3));
        G = 1/sqrt(fctr)*fft2(G);
        res = G-x;
    end
%     toc       
 case {'image'}
    res = zeros(size(x));
    if a.adjoint
        for n=1:nCoils
            tmpk = squeeze(conj(a.KERNEL(:,:,n,:)));
            res(:,:,n) = sum(tmpk.*x,3); 
        end
        res = res - x;
    else
        for n=1:nCoils
            tmpk = a.KERNEL(:,:,:,n);
            res(:,:,n) = sum(tmpk.*x,3); 
        end
        res = res-x;
    end    
end
