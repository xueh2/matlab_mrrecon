classdef SPIRiT_KernelIm
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
        function  res = SPIRiT_KernelIm(kernelIm,method)
        %res = res = SPIRiT_KernelIm(kernelIm,method)
        %	Implementation of the SPIRiT kernel operator
        
            res.KERNEL = kernelIm;
            res.CONJKERNEL = conj(permute(res.KERNEL, [1 2 4 3]));

            % for methods 'fft' and 'image' precompute the image domain kernel by
            % zero-padding and inverse fft

            hasGPU = 1;
            try
                gpuDeviceCount;
            catch
                hasGPU = 0;
            end

            if ( ~hasGPU & (strcmp(method,'fft_GPU')==1 | strcmp(method,'fft_AllGPU')==1) )
                method = 'fft';
            end

            res.adjoint = 0;

            if ( strcmp(method,'fft_GPU')==1 | strcmp(method,'fft_AllGPU')==1 )
                res.KERNELGPU = gpuArray(single(res.KERNEL));
                res.CONJKERNELGPU = gpuArray(single(res.CONJKERNEL));
            else
                res.KERNELGPU = [];
                res.CONJKERNELGPU = [];
            end   

            res.method = method;
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
