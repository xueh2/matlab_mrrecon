classdef SPIRiTMultiOutput
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
        
        KERNEL_I
        CONJKERNEL_I
        identityKernelIm
    end 

    methods 
        function  res = SPIRiTMultiOutput(kernel,method,imSize)
        %res = SPIRiTMultiOutput(kernel [,method, imSize])
        %	Implementation of the SPIRiT kernel operator for multi-output
        %
        %   
        %   Constructor inputs:
        %       kernel: [k_fe, k_pe, srcCha, dstCha, o_fe, o_pe] the spirit 2D convolution kernel.
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

            if nargin < 2
                method = 'fft';
            end

            if strcmp(method,'fft')==1 & nargin < 3
                error('must provide image size');
            end

            % if the multiple output kernel is used, then the input kernel needs expansion
            k_fe = size(kernel, 1);
            k_pe = size(kernel, 2);
            srcCha = size(kernel, 3);
            dstCha = size(kernel, 4);
            o_fe = 1;
            o_pe = 1;

            if ( length(size(kernel)) > 4 )
                o_fe = size(kernel, 5);
                o_pe = size(kernel, 6);
            end

            if ( mod(k_fe,2) ~= 1 )
                error('mod(k_fe,2) ~= 1');
            end

            if ( mod(k_pe,2) ~= 1 )
                error('mod(k_pe,2) ~= 1');
            end

            if ( mod(o_fe,2) ~= 1 )
                error('mod(o_fe,2) ~= 1');
            end

            if ( mod(o_pe,2) ~= 1 )
                error('mod(o_pe,2) ~= 1');
            end

            if ( o_fe > k_fe )
                error('o_fe > k_fe');
            end

            if ( o_pe > k_pe )
                error('o_pe > k_pe');
            end
            
%             newKernel = zeros([2*k_fe-1 2*k_pe-1 srcCha dstCha o_fe o_pe]);
%             
%             hKfe = floor(k_fe/2);
%             hKpe = floor(k_pe/2);
% 
%             hOfe = floor(o_fe/2);
%             hOpe = floor(o_pe/2);
% 
%             for ope = -hOpe:hOpe
%                 for ofe = -hOfe:hOfe
%                     for pe = -hKpe:hKpe
%                         for fe = -hKfe:hKfe
%                             ife = fe-ofe+k_fe;
%                             ipe = pe-ope+k_pe;
%                             newKernel(ife, ipe, :,:,ofe+hOfe+1,ope+hOpe+1) = kernel(fe+hKfe+1,pe+hKpe+1,:,:,ofe+hOfe+1,ope+hOpe+1);                            
%                         end
%                     end                    
%                 end
%             end
%             
%             KERNEL_Im = zeros([imSize(1) imSize(2) srcCha dstCha o_fe o_pe]);
%             for ope = -hOpe:hOpe
%                 for ofe = -hOfe:hOfe            
%                     for n=1:dstCha
%                         aKernel = newKernel(end:-1:1, end:-1:1,:,n, ofe+hOfe+1, ope+hOpe+1);
%                         KERNEL_Im(:,:,:,n,ofe+hOfe+1,ope+hOpe+1) = ifft2c( zpad( aKernel*sqrt(imSize(1)*imSize(2)), imSize(1),imSize(2),srcCha ) );
%                     end
%                 end
%             end
%             
%             KERNEL = zeros([imSize(1) imSize(2) srcCha dstCha]);
%             for ope = 1:o_pe
%                 for ofe = 1:o_fe
%                     KERNEL(:,:,:,:) = KERNEL(:,:,:,:) + KERNEL_Im(:,:,:,:, ofe, ope);
%                 end
%             end
%             KERNEL = KERNEL ./ (o_fe*o_pe);
%             CONJKERNEL = conj(permute(KERNEL, [1 2 4 3]));

            res.KERNEL = convertSPIRiTKSpaceKernel2ImageSpace(kernel, imSize);
            res.CONJKERNEL = conj(permute(res.KERNEL, [1 2 4 3]));
            
            kfe = size(kernel, 1);
            kpe = size(kernel, 2);
            cha = size(kernel, 3);
            
            hkfe = ceil(kfe/2);
            hkpe = ceil(kpe/2);
                       
            %% identity kernel
            identityKernel = zeros(kfe,kpe,cha,cha);
            for d=1:cha
                identityKernel(hkfe, hkpe, d, d) = 1.0;
            end
            res.identityKernelIm = convertSPIRiTKSpaceKernel2ImageSpace(identityKernel, imSize);
           
            res.KERNEL_I = res.KERNEL;
            % res.KERNEL_I = res.KERNEL_I - res.identityKernelIm;            
            res.CONJKERNEL_I = conj(permute(res.KERNEL_I, [1 2 4 3]));
            
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
