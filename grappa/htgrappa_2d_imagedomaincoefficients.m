
function CoeffInIm = htgrappa_2d_imagedomaincoefficients(grappaCoeff, header, kernalLengthFE, blockCoord, missingLines)
%
% CoeffinIm = htgrappa_2d_imagedomaincoefficients(grappaCoeff, header, kernalLengthFE)
%
% ------------------------ input ------------------------
% grappaCoeff : grappa kernel in kspace, [NumberOfWeightsPerMissingLinePerCoil * ( (header.subsampling_factor-1)*header.num_coils)]
% header: information header
% header.blocks is the block size for grappaCoeff
% header.Nfe : number of frequency encoding points
% header.Npe : number of full phase encoding lines
% kernalLengthFE : for 2D grappa kernel, the number of kernel size (number of points used) along the frequency encoding direction; normall (3,5,7 etc)
% blockCoord : define the block structure used by the grappa kernel
% assume the missed k-space lines in a block have the line numbers [1:(header.subsampling_factor-1)], blockCoord gives the line numbers of all lines in this block
% e.g. for reduction factor 4, block size 4, blockCoord can be [-4 0 4 8]
%
% E.g., if we have 12 coil, reduction factor is 4, header.blocks is 4, kernalLengthFE is 5, 
% then grappaCoeff should have the size 240*36
% the 36 columns store the grappa kernels of every missing lines of every coil 
% first column : 1st missing line for coil 1; 2nd : 2nd missing line for coil 1; 3rd : 3rd missing line for coil 1; 4th : 1st missing line for coil 2 etc.
% for every column, the 240 (20 weights per coil * 12coils) weights are stored as 4*12*5 arrays: 
% 1st dimension is line index, corresponding to block lines 1:4
% 2nd dimension is coil, corresponding to 12 coils 1:12
% 3rd dimension is FE shift, corresponds to FE offsets -2:2
% 
% ------------------------ output ------------------------
% CoeffinIm : grappa kernel in image domain [lenFE*lenPE*header.num_coils*header.num_coils], -N/2 to N/2
% e.g. CoeffinIm(:, :, 3, 2) is the 3rd image domain grappa kernel for coil 2
% for every coil, header.num_coils image domain grappa kernel is generated
% 
% Hui Xue
% Jan 24, 2011
%
% References:   HTGRAPPA: real-time B1-weighted image domain TGRAPPA reconstruction MRM 61(6) 2009
% ========================================================

numOfCoils = header.num_coils;
reductionFactor = header.subsampling_factor;
blocks = header.blocks;
halfBlock = floor(kernalLengthFE/2);
Nfe = header.Nfe;
Npe = header.Npe;
lenOfGrappaKernel = blocks * kernalLengthFE;

s = size(grappaCoeff);
if ( s(1) ~= numOfCoils*lenOfGrappaKernel )
    'Error ! the size of grappa kernel array is no correct';
end

if ( s(2) ~= numOfCoils*(reductionFactor-1) )
    'Error ! the size of grappa kernel array is no correct';
end

CoeffInIm = zeros([Nfe Npe numOfCoils numOfCoils]);

% adjust the order of grappa weights to be FE*PE*Coil by Coil*numOfMissingLines
grappaCoeff = adjustWeightsOrder(grappaCoeff, kernalLengthFE, blocks, numOfCoils);

for c = 1:numOfCoils % for every coil
    
    % kernel for coil c
    currKernel = grappaCoeff(:, (c-1)*(reductionFactor-1)+1:c*(reductionFactor-1));    
    
    for im = 1:numOfCoils % for every weight image of this coil (a total of numOfCoils weight images)
    
        % first get the weights for each missing line
        weightsForMissingLines = currKernel( (im-1)*lenOfGrappaKernel+1:im*lenOfGrappaKernel, : );
        
        % fill the conv kernels
        [convKernel, kfe, kpe, feIndexes, peIndexes] = computeConvKernelIndexes(weightsForMissingLines, halfBlock, blockCoord, missingLines, reductionFactor);
                        
        if ( im ~= c)
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
        CoeffInIm(:, :, im, c) = ifftshift(ifft2(CoeffInIm(:, :, im, c)));
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
                        convKernel(dx+halfBlock+1, dy-ymin+1) = weightsForMissingLines(j*feKernelLen+i+1, r);
                    end
                end
            end
        end
        
        % the central point of convKernel is 1
        convKernel(halfBlock+1, -ymin+1) = 1.0;        
    end
    % ------------------------------------------------------------------------------------
    
    % ------------------------------------------------------------------------------------
    % subroutine, adjust the order of grappa weights
    % ------------------------------------------------------------------------------------
    function grappaKernelFePeCoil = adjustWeightsOrder(grappaKernel, lenFE, lenPE, lenCoil)
        % grappaKernel is PE*Coil*PE by Coil*numOfMissingLines array
        % grappaKernelFePeCoil are FE*PE*Coil by Coil*numOfMissingLines array
        s = size(grappaKernel);
        grappaKernel2 = reshape(grappaKernel, [lenPE lenCoil lenFE s(2)]);
        grappaKernel2 = permute(grappaKernel2, [3 1 2 4]);
        grappaKernelFePeCoil = reshape(grappaKernel2, [s(1) s(2)]);        
    end
end
