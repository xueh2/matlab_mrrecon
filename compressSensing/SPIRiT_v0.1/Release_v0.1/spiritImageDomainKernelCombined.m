
function kernelCombined = spiritImageDomainKernelCombined(KERNELS2D, KERNELD2S)
% compute the combined kernel KERNELS2D*KERNELD2S

srcCha = size(KERNELS2D, 3);
dstCha = size(KERNELS2D, 4);
KERNEL = zeros(size(KERNELS2D, 1), size(KERNELS2D, 2), dstCha, dstCha);
for d=1:dstCha
    for dprime=1:dstCha                    
        for s=1:srcCha
            KERNEL(:,:,d,dprime) = KERNEL(:,:,d,dprime) + KERNELD2S(:,:,d,s).*KERNELS2D(:,:,s,dprime);
        end                    
    end
end

kernelCombined = KERNEL;