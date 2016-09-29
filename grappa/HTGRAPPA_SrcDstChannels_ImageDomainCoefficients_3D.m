
function CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients_3D(Coef, option, Nfe, Npe, Npar, reductionFactor, reductionFactorPAR)
%
% CoeffInIm = HTGRAPPA_SrcDstChannels_ImageDomainCoefficients_3D(Coef, header, kernalLengthFE, blockCoord, missingLines)
%
% ------------------------ input ------------------------
% Coef : grappa kernel in kspace, computed using Coef = GRAPPA_SrcDstChannels_Kernel_3D(kspace, acs, option, thresReg)
%
% option.KernelSize, [k_fe k_pe k_par], kernel size along fe and pe and par direction 
% option.KernelPattern, kernel lines along PE, e.g. for R4, it could be [-4 0 4 8]
% option.KernelPatternPAR, kernel lines along PAR, e.g. for R4, it could be [-4 0 4 8]
% option.OutPattern, indexes of lines to be generated, e.g. for R4, it could be [1 2 3] 
% option.OutPatternPAR, indexes of lines to be generated, e.g. for R4, it could be [1 2 3] 
% 
% the size of Coef is [KernelSize(1) KernelSize(2) KernelSize(3) srcCha dstCha length(OutPattern) length(OutPatternPAR)]
% the order of stored coefficients are [fe pe par srcCha dstCha ope opar]
%
% ------------------------ output ------------------------
% CoeffinIm : grappa kernel in image domain [Nfe*Npe*Npar*srcCha*dstCha]
% e.g. CoeffinIm(:, :, :, 3, 2) is the 3rd image domain grappa kernel for destination channel 2
% for every coil, srcCha image domain grappa kernel is generated
% 
% Hui Xue
% Aug 7, 2013
%
% References:   HTGRAPPA: real-time B1-weighted image domain TGRAPPA reconstruction MRM 61(6) 2009
% ========================================================

s = size(Coef);
srcCha = s(4);
dstCha = s(5);

blockCoord = option.KernelPattern;
blockCoordPAR = option.KernelPatternPAR;
halfBlock = floor(option.KernelSize(1)/2);

fitAcquiredLine = false;
if ( ~isempty(find(option.OutPattern==0)) & ~isempty(find(option.OutPatternPAR==0)) )
    fitAcquiredLine = true;
end

CoeffInIm = zeros([Nfe Npe Npar srcCha dstCha], 'single');
CoeffInIm = complex(CoeffInIm, CoeffInIm);

% fill the conv kernels
[convKernel, kfe, kpe, kpar, feIndexes, peIndexes, parIndexes] = computeConvKernelIndexes3D(Coef, halfBlock, blockCoord, blockCoordPAR, option.OutPattern, reductionFactor, option.OutPatternPAR, reductionFactorPAR);
        
if ( ~fitAcquiredLine )
    for c = 1:dstCha % for every destination coil    
        for im = 1:srcCha % for every src channel
            if ( im~=c )
                indfe = find(feIndexes==0);
                indpe = find(peIndexes==0);
                indpar = find(parIndexes==0);
                convKernel(indfe, indpe, indpar, im, c) = 0;
            end
        end
    end
end

% fill the CoeffInIm using the convKernel
% perform zero-filling
for par = parIndexes(1):parIndexes(end)
    for pe = peIndexes(1):peIndexes(end)
        for fe = feIndexes(1):feIndexes(end)                
            x = mod(fe, Nfe);
            y = mod(pe, Npe);
            z = mod(par, Npar);
            CoeffInIm(x+1, y+1, z+1, :, :) = convKernel(fe-feIndexes(1)+1, pe-peIndexes(1)+1, par-parIndexes(1)+1, :, :);
        end
    end
end
        
% get the image domain grappa kernel
for c = 1:dstCha % for every destination coil    
    for im = 1:srcCha % for every src channel                            
        CoeffInIm(:, :, :, im, c) = Nfe*Npe*Npar*ifftshift(ifftn(CoeffInIm(:, :, :, im, c)));
    end
end

    % ------------------------------------------------------------------------------------
    % subroutine, compute the convolution kernel and corresponding k-space indexes
    % ------------------------------------------------------------------------------------
    function [convKernel, kx, ky, kz, xIndexes, yIndexes, zIndexes] = computeConvKernelIndexes3D(weightsForMissingLines, halfBlock, blockCoord, blockCoordPAR, missingLines, reductionFactor, missingLinesPAR, reductionFactorPAR)
        % kx: frequency encoding
        % ky: phase encoding       
        % kz: partition encoding
        
        ymin = min(missingLines) - max(blockCoord);
        ymax = max(missingLines) - min(blockCoord);

        zmin = min(missingLinesPAR) - max(blockCoordPAR);
        zmax = max(missingLinesPAR) - min(blockCoordPAR);

        xIndexes = -halfBlock:halfBlock;
        yIndexes = ymin:ymax;
        zIndexes = zmin:zmax;
        
        feKernelLen = length(xIndexes);
        
        totalNumOfConvKernelPoints = feKernelLen * length(blockCoord) * reductionFactor * length(blockCoordPAR) * reductionFactorPAR;
        
        kx = zeros([1 totalNumOfConvKernelPoints]);
        ky = zeros([1 totalNumOfConvKernelPoints]);
        kz = zeros([1 totalNumOfConvKernelPoints]);
        
        sCha = size(weightsForMissingLines, 4);
        dCha = size(weightsForMissingLines, 5);
        
        % note, the unit of PE is deltaKy
        convKernel = zeros([length(xIndexes) length(ymin:1:ymax) length(zmin:1:zmax) sCha dCha]);

        % if the acquired lines will not be fitted
        if ( (numel(missingLines)<reductionFactor) & (numel(missingLinesPAR)<reductionFactorPAR) )
            ind = 1;
            for x=-halfBlock:halfBlock
                dx = 0 - x;
                                        
                for y = blockCoord(1):reductionFactor:blockCoord(end)
                    for r = 1:reductionFactor
                        if ( r < reductionFactor ) 
                            dy = missingLines(r) - y;
                        else
                            dy = missingLines(end)+1 - y;
                        end
                        
                        for z = blockCoordPAR(1):reductionFactorPAR:blockCoordPAR(end)
                            for rz = 1:reductionFactorPAR

                                if ( rz < reductionFactorPAR ) 
                                    dz = missingLinesPAR(rz) - z;
                                else
                                    dz = missingLinesPAR(end)+1 - z;
                                end
                                
                                kx(ind) = dx;
                                ky(ind) = dy;
                                kz(ind) = dz;
                                ind = ind + 1;

                                i = x + halfBlock ;
                                j = (y - blockCoord(1))/reductionFactor;
                                k = (z - blockCoordPAR(1))/reductionFactorPAR;
                                if ( (r<reductionFactor) | (rz<reductionFactorPAR) )
                                    convKernel(dx+halfBlock+1, dy-ymin+1, dz-zmin+1, :, :) = weightsForMissingLines(i+1, j, k, :, :, r, rz);
                                end
                            end
                        end
                    end
                end
            end
            % the central point of convKernel is 1 unless the acquired lines will be fitted as well
            convKernel(halfBlock+1, -ymin+1, -zmin+1) = 1.0;
        else
            ind = 1;
            for x=-halfBlock:halfBlock
                dx = 0 - x;
                
                for y = blockCoord(1):reductionFactor:blockCoord(end)
                    for r = 1:reductionFactor % go through every missing lines

                        dy = missingLines(r) - y;
                        
                        for z = blockCoordPAR(1):reductionFactorPAR:blockCoordPAR(end)
                            for rz = 1:reductionFactorPAR

                                dz = missingLinesPAR(rz) - z;
                        
                                kx(ind) = dx;
                                ky(ind) = dy;
                                kz(ind) = dz;
                                ind = ind + 1;

                                i = x + halfBlock ;
                                j = (y - blockCoord(1))/reductionFactor;
                                k = (z - blockCoordPAR(1))/reductionFactorPAR;
                                convKernel(dx+halfBlock+1, dy-ymin+1, dz-zmin+1, :, :) = weightsForMissingLines(i+1, j+1, k+1, :, :, r, rz);
                            end
                        end
                    end
                end
            end
        end
    end
    % ------------------------------------------------------------------------------------
end