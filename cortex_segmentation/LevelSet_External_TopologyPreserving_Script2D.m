
clear
cd J:\Neonatal_brain_data\test_TopologyPreserving
%--------------------------------------------------------------------------
% filename = 'testImage.bmp';

filename = 'singlecircle2.bmp';

[x, map] = imread(filename, 'bmp');
x = double(x);
x=~x;

[ysize, xsize] = size(x);
xvoxelsize = 1;
yvoxelsize = 1;

header.xsize = xsize;
header.ysize = ysize;
header.zsize = 1;
header.xvoxelsize = 1;
header.yvoxelsize = 1;
header.zvoxelsize = 1;
header.bytes = 2;

SaveAnalyze(uint32(x), header, 'x.hdr', 'Grey'); 

SignedPressureForce = 2*x-1;
imview(SignedPressureForce, []);

% perform the level set propagation
%---------------------------------------------------------------------------
normalSpeed = SignedPressureForce;
bValue = 0.02;

accuracy = 'medium';
tMax_ReIntialize = 10;
errorMax = 0.05;

resultDir = 'levelset_topologypreserving';

tMax = 35;                    % End time.
% plotSteps = 21;              % How many intermediate plots to produce?
tPlot = 0.2;
factorCFL = 0.5;
reInitialStep = 4;

flag_outside = 0  % target surface must be outside the data0 surface
flag_inside = 0;  % target surface must be inside the data0 surface

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 1;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 0;

displayType = 'contour';
%--------------------------------------------------------------------------
% Integration parameters.
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
% tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;
%--------------------------------------------------------------------------
% Create the grid.
g.dim = 2;
% ==========================
g.min = [-yvoxelsize * (ysize-1)/2.0; -xvoxelsize * (xsize-1)/2.0];
g.max = [yvoxelsize * (ysize-1)/2.0; xvoxelsize * (xsize-1)/2.0];

g.dx = [yvoxelsize; xvoxelsize];

g.bdry = @addGhostExtrapolate;

g = processGrid(g, x);

%---------------------------------------------------------------------------
% create SDF
% center = [45;45];
% center = center + g.min;
% radius = [30];
% SDF = shapeSphere(g, center, radius)
% data0 = SDF;

% center = [66;131];
% center = center + g.min;
% radius = [40];
% SDF1 = shapeSphere(g, center, radius);
% 
% center = [66;37];
% center = center + g.min;
% radius = [30];
% SDF2 = shapeSphere(g, center, radius);
% data0 = shapeunion(SDF1, SDF2);

filename = 'sdf4.bmp';
[x, map] = imread(filename, 'bmp');
x = ~x;
x = double(x);
SDF = CreateApproximated_SDF_2D(x);
data0 = SDF;

SaveAnalyze(data0, header, 'data0.hdr', 'Real'); 

%---------------------------------------------------------------------------

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', factorCFL, 'stats', 'on');

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

%---------------------------------------------------------------------------
% Set up motion in the normal direction.
normalFunc = @termNormal;
normalData.grid = g;
normalData.speed = normalSpeed;
normalData.derivFunc = derivFunc;

%---------------------------------------------------------------------------
curvatureFunc = @termCurvature;
curvatureData.grid = g;
curvatureData.curvatureFunc = @curvatureSecond;
curvatureData.b = bValue;

%--------------------------------------------------------------------------

% Convergence criteria
deltaMax = errorMax * max(g.dx) * prod(g.N);

%---------------------------------------------------------------------------
% Combine components of motion.
if ( curvatureData.b ~= 0 )
    schemeFunc = @termSum;
    schemeData.innerFunc = { normalFunc; curvatureFunc };
    schemeData.innerData = { normalData; curvatureData };
else
    schemeFunc = @termSum;
    schemeData.innerFunc = { normalFunc };
    schemeData.innerData = { normalData };
end
schemeData.data0 = data0;
schemeData.y0 = data0(:);
schemeData.shape = g.shape;

% for max/min operation
schemeData.flag_outside = flag_outside;
schemeData.flag_inside = flag_inside;

% for the topology perserving level set
schemeData.B = zeros(size(data0), 'uint8');
schemeData.LastY = schemeData.y0;
schemeData.LastData = data0;
schemeData.connectivityObject = 4;
schemeData.connectivityBackground = 8;
schemeData = processTopology(schemeData);
%---------------------------------------------------------------------------
% Let the integrator know what function to call.
integratorOptions = odeCFLset(integratorOptions, ...
                              'postTimestep', @CheckTopology2D);
% integratorOptions = odeCFLset(integratorOptions, ...
%                               'postTimestep', @GetIntermediateDigitalTopology);
%--------------------------------------------------------------------------
% Initialize Display
f = figure;

% Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

h = visualizeLevelSet(g, data0, displayType, level, [ 't = ' num2str(t0) ]);

% hold on;
if(g.dim > 1)
  axis(g.axis);
  daspect([ 1 1 1 ]);
end

%--------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
data = data0;
tLastReInitial = t0;

global ntimes
ntimes = 1;

while(tMax - tNow > small * tMax)
    disp(['tNow = ' num2str(tNow) ' ... ']);
    % Reshape data array into column vector for ode solver call.
    y0 = data(:);

    % How far to step?
    tSpan = [ tNow, min(tMax, tNow + tPlot) ];

    % Take a timestep.
    [ t, y, schemeData ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
    tNow = t(end);

    % Get back the correctly shaped data array
    data = reshape(y, g.shape);

%     if ( flag_outside )
%         data = min(data, data0);
%     end
%     
%     if ( flag_inside )
%         data = max(data, data0);
%     end
    % Check for convergence (except for the first loop).
    
%     if( norm(y - y0, 1) < deltaMax )
%         disp('the update of phi is too small ...');
%         break;
%     end
    if(pauseAfterPlot)
        % Wait for last plot to be digested.
        pause;
    end

    % Get correct figure, and remember its current view.
    %figure(f);
    figure(f);
    figureView = view;
    clf;

    % Delete last visualization if necessary.
%     if(deleteLastPlot)
%         delete(h);
%     end

    % Move to next subplot if necessary.
    if(useSubplots)
        plotNum = plotNum + 1;
        subplot(rows, cols, plotNum);
    end

    % Create new visualization.
    h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
%     imtool(data);
    % Restore view.
%     view(figureView);

%     filename = ['data' num2str(tNow) '.hdr'];
%     SaveAnalyze(data, header, filename, 'Real'); 
    
    if ( (tNow - tLastReInitial) >= reInitialStep );
        disp(['tNow = ' num2str(tNow) ': Reinitializing ... ']);
        tLastReInitial = tNow;
        
        data = signedDistanceIterative(g, data, accuracy, tMax_ReIntialize, errorMax);

    end
end

endTime = cputime;
fprintf('Total execution time %g seconds', endTime - startTime);

header.xsize = xsize;
header.ysize = ysize;
header.zsize = 1;
header.xvoxelsize = 1;
header.yvoxelsize = 1;
header.zvoxelsize = 1;
header.bytes = 2;

SaveAnalyze(data, header, 'levelset.hdr', 'Real');