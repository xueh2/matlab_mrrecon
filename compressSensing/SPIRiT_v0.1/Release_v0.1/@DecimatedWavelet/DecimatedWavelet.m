
classdef DecimatedWavelet
    % Class help goes here
    properties
        adjoint
        N
        wname
        ExtM
    end 

    methods 
        function res = DecimatedWavelet(N, wname, ExtM)
            % res = DecimatedWavelet(N, wname, ExtM)
            %
            % implements a wavelet operator. 
            % This is a decimated wavelet operator

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
        mag = coeffNorm(coeff)
    end
end