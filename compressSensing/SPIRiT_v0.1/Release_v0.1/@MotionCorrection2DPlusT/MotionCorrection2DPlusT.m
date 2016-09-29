
% Hui Xue 2012.03.03

classdef MotionCorrection2DPlusT
    % Class help goes here
    properties
        adjoint
        header
        interpolator
        keyFrame
        dx
        dy
        dxInv
        dyInv
    end 

    methods 
        function res = MotionCorrection2DPlusT(interpolator, keyFrame, header, dx, dy, dxInv, dyInv)
            % res = MotionCorrection2DPlusT(N)
            % implements the motion correction operator on 2D+T complex Im
            % interpolator : 'NN', 'BSpline', 'Linear'
            res.adjoint = 0;
            res.header = header;            
            res.interpolator = interpolator;
            res.keyFrame = keyFrame;
            
            if ( nargin < 7 )
                res.dx = [];
                res.dy = [];
                res.dxInv = [];
                res.dyInv = [];                
            else
                res.dx = dx;
                res.dy = dy;
                res.dxInv = dxInv;
                res.dyInv = dyInv;
            end
        end

        res = mtimes(a,b)
        res = times(a,b)
        res = ctranspose(a)
        [res, a] = performMoCo(a, im, params)
        a = computeConcatenatedDeformationField(a, dx, dy, dxInv, dyInv)
    end
end