
% Hui Xue 2012.02.20

classdef RedundantHarrDWT3D
    % Class help goes here
    properties
        adjoint
        N
        temporalSoftThreshScaleRatio
    end 

    methods 
        function res = RedundantHarrDWT3D(N)
            % res = RedundantHarrDWT3D(N)
            %
            % implements the redundant harr wavelet operator using MrFtk
            % This is an undecimated wavelet operator

            res.adjoint = 0;
            res.N = N;
            res.temporalSoftThreshScaleRatio = 1;
        end

        res = mtimes(a,b)
        res = times(a,b)
        res = ctranspose(a)
        res = softThres(y, t)
        coeff = scaleTemporalCoeff(a,coeff,temporalScalingFactor)
    end
    
    methods(Static)        
        [coeffNorm, totalNorm] = coeffNorm(coeff)
        coeff = divideByNorm(coeff, coeffNorm, p, mu)
    end
end