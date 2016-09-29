function kernelIm = convertSPIRiTKSpaceKernel2ImageSpace(kernel, imSize)
% This method conver the kspace spirit kernel to image space
% kernel: [k_fe, k_pe, srcCha, dstCha, o_fe, o_pe] the spirit 2D convolution kernel with multiple outputs
% imSize : [Nfe Npe], the size of image space kernel
% kernelIm : image space kernels

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

newKernel = zeros([2*k_fe-1 2*k_pe-1 srcCha dstCha o_fe o_pe]);

hKfe = floor(k_fe/2);
hKpe = floor(k_pe/2);

hOfe = floor(o_fe/2);
hOpe = floor(o_pe/2);

for ope = -hOpe:hOpe
    for ofe = -hOfe:hOfe
        for pe = -hKpe:hKpe
            for fe = -hKfe:hKfe
                ife = fe-ofe+k_fe;
                ipe = pe-ope+k_pe;
%                 kernel(fe+hKfe+1,pe+hKpe+1,1,1,ofe+hOfe+1,ope+hOpe+1)
                
                newKernel(ife, ipe, :,:,ofe+hOfe+1,ope+hOpe+1) = kernel(fe+hKfe+1,pe+hKpe+1,:,:,ofe+hOfe+1,ope+hOpe+1);                            
            end
        end                    
    end
end

% KERNEL_Im = zeros([imSize(1) imSize(2) srcCha dstCha o_fe o_pe]);
% for ope = -hOpe:hOpe
%     for ofe = -hOfe:hOfe            
%         for n=1:dstCha
%             aKernel = newKernel(end:-1:1, end:-1:1,:,n, ofe+hOfe+1, ope+hOpe+1);
%             KERNEL_Im(:,:,:,n,ofe+hOfe+1,ope+hOpe+1) = zpad(aKernel*sqrt(imSize(1)*imSize(2)), imSize(1),imSize(2),srcCha );
%         end
%     end
% end
% 
% KERNEL_Im = sum(KERNEL_Im, 6);
% KERNEL_Im = sum(KERNEL_Im, 5);
% fctr = imSize(1)*imSize(2);
% kernelIm = sqrt(fctr)*fftshift(ifft2(ifftshift(KERNEL_Im))) ./ (o_fe*o_pe);

newKernel = sum(newKernel, 6);
newKernel = sum(newKernel, 5);

KERNEL_Im = zeros([imSize(1) imSize(2) srcCha dstCha]);            
for n=1:dstCha
    aKernel = newKernel(end:-1:1, end:-1:1,:,n);
    KERNEL_Im(:,:,:,n) = zpad(aKernel*sqrt(imSize(1)*imSize(2)), imSize(1),imSize(2),srcCha );
end

fctr = imSize(1)*imSize(2);
kernelIm = sqrt(fctr)*fftshift(ifft2(ifftshift(KERNEL_Im))) ./ (o_fe*o_pe);
