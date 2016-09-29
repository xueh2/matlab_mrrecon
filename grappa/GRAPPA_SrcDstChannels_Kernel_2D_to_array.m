

function CoefArray = GRAPPA_SrcDstChannels_Kernel_2D_to_array(Coef, kfe, kpe, srcCha, dstCha, ope)
% function Coef = GRAPPA_SrcDstChannels_Kernel_2D_to_array(Coef, CoefArray)
% Coef : grappa kernel estimated
% the size of Coef is [KernelSize(1)*KernelSize(2)*srcCha length(OutPattern)*dstCha]
% e.g., for a [5 4] kernel, 32 src channels and 4 dst channels, R4, missing lines [1 2 3]
% the coefficient size is [5*4*32 3*4]
% the order of stored coefficients are [fe pe srcCha missingLine dstCha]
%
% CoefArray: an array of [fe pe srcCha dstCha ope]

CoefArray = zeros(kfe, kpe, srcCha, dstCha, ope);

for o=1:ope
    for d=1:dstCha
        for s=1:srcCha
            for pe=1:kpe
                for fe=1:kfe
                                        
                    CoefArray(fe, pe, s, d, o) = Coef(fe+(pe-1)*kfe+(s-1)*kfe*kpe, o+(d-1)*ope);
                    
                end
            end
        end
    end
end
