
% Hui Xue 2011

classdef PRUNO
    % Implement the PRUNO MR recon
    % ref:Parallel reconstruction using null operations, Jian Zhang, Chunlei Liu, Michael E. Moseley
    
    properties
        adjoint % adjoint operator
        Wd      % width of kernel
        r       % number of kernels selected to build the null space
        kernel  % pruno kernel [row col nCoil r]
        kernelIm % image domain pruno kernel [Nfe Npe nCoil r]
        kernelComposite % composited kernel [Nfe Npe nCoil nCoil]
    end 

    methods 
        function res = PRUNO(Wd, r)
            % res = PRUNO(Wd, r)
            % This PURNO operator implements N'N defined in the paper
            %
            % implements a purno operator. 

            res.adjoint = 0;
            res.Wd = Wd;
            res.r = r;
        end

        res = mtimes(a,b)
        res = times(a,b)
        res = ctranspose(a)
        res = transpose(a)
    end
    
    methods(Static)
        kernel = prunoCalibration(kCalib,kSize,thresEigenValue, r, show)
        kernelIm = convertKernelToImage(kernel, imSize, show, voxelsize, centre, width)
        kernelComposite = computeCompositeNullKernel(kernelIm, show, voxelsize, centre, width)
    end
end