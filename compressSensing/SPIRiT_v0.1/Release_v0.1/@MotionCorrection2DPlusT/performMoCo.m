function [res, a] = performMoCo(a,b,params)

if isa(a,'MotionCorrection2DPlusT') == 0
    error('In performMoCo(a,b,params), a must be MotionCorrection2DPlusT operator');
end

if ( nargin < 3 )
    params.strategy = 'Direct';
    params.inverse = 1; 
    params.initial = 0; 
    params.numOfPre = 0; 
    params.iters = [100 100 100]; 
    params.sigma = 24.0; 
    params.neighbor = 2.0; 
    params.stepDiv =  3.0; 
    params.moreIterInv = 1; 
    params.algo = 'GLCC'; 
    params.volumePreserving = 0;
    % every refineMoCoDuringReconIter iterations, the moco is performed and deformation fields are refined, if -1, no moco refinement will be performed
    params.numOfIterToRefineMoCoDuringRecon = 2; 
    % every refinement for moco, the sigma is decreased by this ratio
    params.sigmaDivRatio = 1.5; 
end

[res, a.dx, a.dy, a.dxInv, a.dyInv] = PerformTemporalMotionCorrectionComplex(b, a.header, a.keyFrame, params.strategy, params.inverse, ...
                        params.initial, params.numOfPre, params.iters, params.sigma, params.neighbor, params.stepDiv, params.moreIterInv, ... 
                        params.algo, params.volumePreserving, a.interpolator);

a.dx = single(a.dx);
a.dy = single(a.dy);
a.dxInv = single(a.dxInv);
a.dyInv = single(a.dyInv);
