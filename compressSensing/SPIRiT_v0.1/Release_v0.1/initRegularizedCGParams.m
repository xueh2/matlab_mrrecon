function res = initRegularizedCGParams()
% param = initRegularizedCGParams()
%
% function returns a structure with the entries that are needed for the reconstruction.

res.xfmWeight = 0.0025;   % transform l1 penalty
res.tikWeight = 0;
res.TVWeight = 1.0e-004;

res.continuationStep = 5;
res.wavWeightRatio = 1;
                    
res.Itnlim = 70;	% default number of iterations
res.gradToll = 1e-30;	% step size tollerance stopping criterea
res.objToll = 1e-4;	% objection function tollerance stopping criterea
res.l1Smooth = 1e-15;	% smoothing parameter of L1 norm
res.pNorm = 1;  % type of norm to use (i.e. L1 L2 etc)

% line search parameters
% res.lineSearchItnlim = 50;
% res.lineSearchAlpha = 0.005;
% res.lineSearchBeta = 0.75;
% res.lineSearchT0 = 1 ; % step size to start with

res.lineSearchItnlim = 150;
res.lineSearchAlpha = 0.005;
res.lineSearchBeta = 0.75;
res.lineSearchT0 = 1 ; % step size to start with

res.show = 0;
% res.continuationRatio = 1;

res.TV = TVOP;


