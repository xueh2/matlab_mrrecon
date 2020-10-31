
function [SDF, header] = SignedDistranceTransform(InputLevelSet_File, OutputSDF_File, accuracy, tMax_ReIntialize, errorMax)
% transform the input level set function to the signed distance function

p = dir(InputLevelSet_File);

if ( isempty(p)==1 )
    error('Can not find InputLevelSet_File ...');
    return;
end

[data, header] = LoadAnalyze(InputLevelSet_File, 'Real');
data = double(data);

if ( exist('accuracy')==0 )
    accuracy = 'high';
end    

if ( exist('tMax_ReIntialize')==0 )
    tMax_ReIntialize = 100;
end    

if ( exist('errorMax')==0 )
    errorMax = 0.05;
end

xsize = header.xsize;
ysize = header.ysize;
zsize = header.zsize;

xvoxelsize = header.xvoxelsize;
yvoxelsize = header.yvoxelsize;
zvoxelsize = header.zvoxelsize;

disp('Signed Distance Function ...');

%---------------------------------------------------------------------------
% Integration parameters.
%--------------------------------------------------------------------------
% Create the grid.
g = [];
g.dim = 3;
% ==========================
% g.min = [-xvoxelsize * (xsize-1)/2.0; -yvoxelsize * (ysize-1)/2.0; -zvoxelsize * (zsize-1)/2.0];
% g.max = [xvoxelsize * (xsize-1)/2.0; yvoxelsize * (ysize-1)/2.0; zvoxelsize * (zsize-1)/2.0];

g.min = [-yvoxelsize * (ysize-1)/2.0; -xvoxelsize * (xsize-1)/2.0; -zvoxelsize * (zsize-1)/2.0];
g.max = [yvoxelsize * (ysize-1)/2.0; xvoxelsize * (xsize-1)/2.0; zvoxelsize * (zsize-1)/2.0];

g.dx = [header.xvoxelsize; header.yvoxelsize; header.zvoxelsize];

g.bdry = @addGhostExtrapolate;

g = processGrid(g, data);

%---------------------------------------------------------------------------
% Choose approximations at appropriate level of accuracy.
%   Same accuracy is used by both components of motion.
switch(accuracy)
 case 'low'
  derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

startTime = cputime;

disp(['tranforming ...']);

SDF = signedDistanceIterative(g, data, accuracy, tMax_ReIntialize, errorMax);

endTime = cputime;
fprintf('Total execution time %g seconds', endTime - startTime);

if ( isempty(OutputSDF_File)==0 )
    SaveAnalyze(single(SDF), header, OutputSDF_File, 'Real');
end

return;
