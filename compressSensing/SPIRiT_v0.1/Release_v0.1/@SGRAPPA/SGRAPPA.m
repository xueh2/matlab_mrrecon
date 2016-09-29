
% Hui Xue, July 2011

classdef SGRAPPA
    % Implement the self-consistent grappa MR recon
    % ref : originally proposed
    
    properties
        adjoint % adjoint operator
        idxacq  % index for availabel kspace points (in this case, original missing points)
        kernel  % self-consistent grappa kernel [row col nCoil nCoil]
        kernelIm % image domain pruno kernel [Nfe Npe nCoil nCoil]
    end 

    methods 
        function res = SGRAPPA(idxacq)
            % res = SGRAPPA()
            % This SGRAPPA operator implements sgrappa
            res.adjoint = 0;
            res.idxacq = idxacq;
        end

        res = mtimes(a,b)
        res = times(a,b)
        res = ctranspose(a)
        res = transpose(a)
    end
    
    methods(Static)
        kernel = sgrappaCalibration(kCalib,kSize,acquiredLines,thresReg)
        kernelIm = sgrappaKernelToImage(kernel, imSize, show, voxelsize, centre, width)
    end
end
