
% Hui Xue 2011

classdef UndecimatedWavelet
    % Class help goes here
    properties
        adjoint
        N
        wname
        ExtM
    end 

    methods 
        function res = UndecimatedWavelet(N, wname, ExtM)
            % res = UndecimatedWavelet(N, wname, ExtM)
            %
            % implements a wavelet operator. 
            % This is an undecimated wavelet operator

            res.adjoint = 0;
            res.N = N;
            res.wname = wname;
            res.ExtM = ExtM;
        end

        res = mtimes(a,b)
        res = times(a,b)        
    end
    
    methods(Static)
        x = softThresh(y,t)
        [coeffNorm, totalNorm] = coeffNorm(coeff)
        coeff = divideByNorm(coeff, coeffNorm, p, mu)
    end
end