
% Hui Xue 2011

classdef RedundantHarrDWT2D
    % Class help goes here
    properties
        adjoint
        N
    end 

    methods 
        function res = RedundantHarrDWT2D(N)
            % res = RedundantHarrDWT2D(N)
            %
            % implements the redundant harr wavelet operator using MrFtk
            % This is an undecimated wavelet operator

            res.adjoint = 0;
            res.N = N;
        end

        res = mtimes(a,b)
        res = times(a,b)
        res = ctranspose(a)
        x = softThresh(a,y,t,gFactor)
    end
    
    methods(Static)        
        [coeffNorm, totalNorm] = coeffNorm(coeff)
        coeff = divideByNorm(coeff, coeffNorm, p, mu)
    end
end