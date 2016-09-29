classdef SPIRiT2DPlusT
    % Class help goes here
    properties
        method;
        adjoint
        numOfRep
        SpiritGOPs;
    end 

    methods 
        function  res = SPIRiT2DPlusT(GOPs)
        %res = SPIRiT2DPlusT(GOPs,method,imSize,numOfRep)
        %	Implementation of the SPIRiT 2D plus T kernel operator
            res.numOfRep = numel(GOPs);            
            res.SpiritGOPs = GOPs;
            res.adjoint = 0;
            res.method = GOPs{1}.method;
        end
        
        res = mtimes(a,b)
        res = times(a,b)
        res = ctranspose(a)
        res = transpose(a)
    end
    
    methods(Static)
        
    end
end
