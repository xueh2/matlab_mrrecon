
function CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients_2(Coef, option, Nfe, Npe, reductionFactor)
%
% CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients(Coef, header, kernalLengthFE, blockCoord, missingLines)
%
% ------------------------ input ------------------------
% Coef : grappa kernel in kspace, computed using Coef = GRAPPA_SrcDstChannels_Kernel_2D(kspace, acs, option, thresReg)
%
% option.KernelSize, [k_fe k_pe], kernel size along fe and pe direction 
% option.KernelPattern, kernel lines along PE, e.g. for R4, it could be [-4 0 4 8]
% option.OutPattern, indexes of lines to be generated, e.g. for R4, it could be [1 2 3] 
% 
% the size of Coef is [KernelSize(1)*KernelSize(2)*srcCha length(OutPattern)*dstCha]
% e.g., for a [5 4] kernel, 32 src channels and 4 dst channels, R4, missing lines [1 2 3]
% the coefficient size is [5*4*32 3*4]
% the order of stored coefficients are [fe pe srcCha missingLine dstCha]
% that is, for Coef(1:20,1), the first 5 elements are first column of kernel for srcCha=1
% for Coef(1:20,1), the next 5 elements are second column of kernel for srcCha=1
% for Coef(21:40,1), the first 5 elements are first column of kernel for srcCha=2
% Coef(:, 1:3) are kernels for first dst channel, Coef(:, 4:6) are kernels for the second dst channel etc.

% ------------------------ output ------------------------
% CoeffinIm : grappa kernel in image domain [Nfe*Npe*srcCha*dstCha]
% e.g. CoeffinIm(:, :, 3, 2) is the 3rd image domain grappa kernel for destination channel 2
% for every coil, srcCha image domain grappa kernel is generated
% 
% Hui Xue
% Oct 10, 2011
%
% References:   HTGRAPPA: real-time B1-weighted image domain TGRAPPA reconstruction MRM 61(6) 2009
% ========================================================

s = size(Coef);
srcCha = size(Coef, 3);
dstCha = size(Coef, 4);

if ( s(1) ~= srcCha*prod(option.KernelSize) )
    'Error ! the size of grappa kernel array is no correct';
end

if ( s(2) ~= dstCha*numel(option.OutPattern) )
    'Error ! the size of grappa kernel array is no correct';
end

blockCoord = option.KernelPattern;
halfBlock = floor(option.KernelSize(1)/2);

numOfMissingLines = numel(option.OutPattern);
lenOfGrappaKernel = option.KernelSize(1)*option.KernelSize(2);
fitAcquiredLine = ~isempty(find(option.OutPattern==0));

CoeffInIm = zeros([Nfe Npe srcCha dstCha]);

for c = 1:dstCha % for every destination coil
    
    % kernel for dst coil c
    currKernel = Coef(:, :, :, c, :);    
    
    for im = 1:srcCha % for every src channel
    
        % first get the weights for each missing line
        weightsForMissingLines = currKernel(:,:,im, 1, :);
        
        % fill the conv kernels
        [convKernel, kfe, kpe, feIndexes, peIndexes] = computeConvKernelIndexes(weightsForMissingLines, halfBlock, blockCoord, option.OutPattern, reductionFactor);
                        
        if ( ~fitAcquiredLine & (im~=c) )
            indfe = find(feIndexes==0);
            indpe = find(peIndexes==0);
            convKernel(indfe, indpe) = 0;
        end
               
        % fill the CoeffInIm using the convKernel
        % CoeffInIm: [0:Nfe-1] * [0:Npe-1]
        % perform zero-filling
        for pe = peIndexes(1):peIndexes(end)
            for fe = feIndexes(1):feIndexes(end)                
                x = mod(fe, Nfe);
                y = mod(pe, Npe);
                CoeffInIm(x+1, y+1, im, c) = convKernel(fe-feIndexes(1)+1, pe-peIndexes(1)+1);
            end
        end
        
        % get the image domain grappa kernel
        CoeffInIm(:, :, im, c) = Nfe*Npe*ifftshift(ifft2(CoeffInIm(:, :, im, c)));
    end
end

    % ------------------------------------------------------------------------------------
    % subroutine, compute the convolution kernel and corresponding k-space indexes
    % ------------------------------------------------------------------------------------
    function [convKernel, kx, ky, xIndexes, yIndexes] = computeConvKernelIndexes(weightsForMissingLines, halfBlock, blockCoord, missingLines, reductionFactor)
        % kx: frequency encoding
        % ky: phase encoding       
               
        ymin = min(missingLines) - max(blockCoord);
        ymax = max(missingLines) - min(blockCoord);

        xIndexes = -halfBlock:halfBlock;
        yIndexes = ymin:ymax;
        
        feKernelLen = length(xIndexes);
        
        totalNumOfConvKernelPoints = feKernelLen * length(blockCoord) * reductionFactor;
        
        kx = zeros([1 totalNumOfConvKernelPoints]);
        ky = zeros([1 totalNumOfConvKernelPoints]);
        
        % note, the unit of PE is deltaKy
        convKernel = zeros([length(xIndexes) length(ymin:1:ymax)]);

        % if the acquired lines will not be fitted
        if ( numel(missingLines) < reductionFactor )
            ind = 1;
            for x=-halfBlock:halfBlock
                for y = blockCoord(1):reductionFactor:blockCoord(end)
                    for r = 1:reductionFactor
                        dx = 0 - x;
                        if ( r < reductionFactor ) 
                            dy = missingLines(r) - y;
                        else
                            dy = missingLines(end)+1 - y;
                        end
                        kx(ind) = dx;
                        ky(ind) = dy;
                        ind = ind + 1;

                        i = x + halfBlock ;
                        j = (y - blockCoord(1))/reductionFactor;
                        if ( r < reductionFactor )
                            convKernel(dx+halfBlock+1, dy-ymin+1) = weightsForMissingLines(i+1, j, 1, 1, r);
                        end
                    end
                end
            end
            % the central point of convKernel is 1 unless the acquired lines will be fitted as well
            convKernel(halfBlock+1, -ymin+1) = 1.0;
        else
            ind = 1;
            for x=-halfBlock:halfBlock
                dx = 0 - x;
                for y = blockCoord(1):reductionFactor:blockCoord(end)
                    for r = 1:reductionFactor % go through every missing lines

                        dy = missingLines(r) - y;
                        
                        kx(ind) = dx;
                        ky(ind) = dy;
                        ind = ind + 1;

                        i = x + halfBlock ;
                        j = (y - blockCoord(1))/reductionFactor;
                        convKernel(dx+halfBlock+1, dy-ymin+1) = weightsForMissingLines(i+1, j+1, 1, 1, r);
                    end
                end
            end
        end
    end
    % ------------------------------------------------------------------------------------
end