classdef SPIRiT3D
    % Class help goes here
    properties
        adjoint
        method
	    KERNEL
        CONJKERNEL        
        kernel
        KERNELGPU
        CONJKERNELGPU
        imSize
    end 

    methods 
        function  res = SPIRiT3D(kernel,method,imSize)
        %res = SPIRiT3D(kernel [,method, imSize])
        %	Implementation of the 3D SPIRiT kernel operator
        %
        %   
        %   Constructor inputs:
        %       kernel: [kx,ky,kz,nCoils,nCoils] the spirit 3D convolution kernel.
        %               See corrMatrix.m, calibrate.m and the demos on how to
        %               generate a kernel.
        %       method: implementation type of the operator. There are three
        %               'conv'  is a k-space convolution implementation (slowest)
        %                       that operates on k-space data and produces k-space data.
        %               'fft'   is an fft based implementation of the k-space
        %                       convolution through image domain multiplication. 
        %                       It operates on k-sapce data and produces k-space data. 
        %               'image' The SPIRiT operator is applied to image space data.
        %                       This is useful when using image-based non-cartesian
        %                       reconstruction. In this case the SPIRiT operator
        %                       operates on image data and produces image data. 
        %
        %       imSize:  Size of the resulting image (only needed for 'fft' and
        %                'image' modes
        %
        %       what is estimated is (G-I) and (G-I)'
        
            if nargin < 2
                method = 'conv';
                KERNEL = [];
                CONJKERNEL = [];
            end

            if strcmp(method,'conv')==1
                KERNEL = [];
                CONJKERNEL = [];
            end


            if strcmp(method,'fft')==1 & nargin < 3
                error('must provide image size');
            end

            % for methods 'fft' and 'image' precompute the image domain kernel by
            % zero-padding and inverse fft

            kernelOnes = kernel;
            kernelOnes(:) = 0;
            
            % only the center is 1
            kfe = size(kernelOnes,1);
            kpe = size(kernelOnes,2);
            kpar = size(kernelOnes,3);
            
            halfKfe = floor(kfe/2)+1;
            halfKpe = floor(kpe/2)+1;
            halfKpar = floor(kpar/2)+1;
            
            srcCha = size(kernelOnes, 4);
            dstCha = size(kernelOnes, 5);
            
            for dst=1:dstCha
                kernelOnes(halfKfe, halfKpe, halfKpar, dst, dst) = 1;
            end                        
            
            kernel = kernel - kernelOnes; % G-I kernel
            
            res.KERNEL = zeros(imSize(1), imSize(2), imSize(3), size(kernel,4), size(kernel,5), 'single');
            res.CONJKERNEL = zeros(imSize(1), imSize(2), imSize(3), size(kernel,5), size(kernel,4), 'single');
            
            if strcmp(method,'fft')==1 | strcmp(method,'image')==1 | strcmp(method,'fft_MrFtk')==1 | strcmp(method,'fft_GPU')==1 | strcmp(method,'fft_AllGPU')==1
                for n=1:size(kernel,5)
                    n
                    res.KERNEL(:,:,:,:,n) = (ifft3c(zpad(kernel(end:-1:1,end:-1:1,end:-1:1,:,n)*sqrt(imSize(1)*imSize(2)*imSize(3)), imSize(1), imSize(2), imSize(3), size(kernel,4))));
                end        
                % res.CONJKERNEL = conj(permute(res.KERNEL, [1 2 3 5 4]));
                for dstCha = 1:size(kernel,5)
                    for srcCha = 1:size(kernel,4)
                        res.CONJKERNEL(:,:,:,dstCha,srcCha) = conj(res.KERNEL(:,:,:,srcCha,dstCha));
                    end
                end
            end

            if strcmp(method,'fft')==0 & strcmp(method,'conv')==0 & strcmp(method,'image')==0 & strcmp(method,'fft_MrFtk')==0 & strcmp(method,'fft_GPU')==0 & strcmp(method,'fft_AllGPU')==0
                eror('no such method');
            end

            hasGPU = 1;
            try
                gpuDeviceCount;
            catch
                hasGPU = 0;
            end

            if ( ~hasGPU & (strcmp(method,'fft_GPU')==1 | strcmp(method,'fft_AllGPU')==1) )
                method = 'fft';
            end

            res.kernel = kernel;
            res.adjoint = 0;

            if ( strcmp(method,'fft_GPU')==1 | strcmp(method,'fft_AllGPU')==1 )
                res.KERNELGPU = gpuArray(single(res.KERNEL));
                res.CONJKERNELGPU = gpuArray(single(res.CONJKERNEL));
            else
                res.KERNELGPU = [];
                res.CONJKERNELGPU = [];
            end   

            res.method = method;
            res.imSize = imSize;
        end
        
        res = mtimes(a,b)
        res = times(a,b)
        res = ctranspose(a)
        res = transpose(a)
        [K Kconv] = getK(G)
        kernel = getKernel(GOP)
        method = getMethod(GOP)
    end
    
    methods(Static)
        
    end
end
