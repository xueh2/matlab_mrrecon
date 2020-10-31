
clear
cd J:\Neonatal_brain_data\test_TopologyPreserving
%--------------------------------------------------------------------------
% signed force
% signed distance function

% singlecircle

xsize = 64;
ysize = 64;
zsize = 40;

xvoxelsize = 1;
yvoxelsize = 1;
zvoxelsize = 1;

header.xsize = xsize;
header.ysize = ysize;
header.zsize = zsize;
header.xvoxelsize = xvoxelsize;
header.yvoxelsize = yvoxelsize;
header.zvoxelsize = zvoxelsize;
header.bytes = 2;

center_binary = [40;40;20];
radius_binary = [16];

center_sdf = [22;22;20];
radius_sdf = [12];

schemeData.connectivityObject = 18;
schemeData.connectivityBackground = 6;

resultDir = 'levelset_topologypreserving_3D_18_6';

%--------------------------------------------------------------------------

accuracy = 'low';
tMax_ReIntialize = 10;
errorMax = 0.05;

mkdir(resultDir);

tMax = 100;                    % End time.
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

displayType = 'surface';
%--------------------------------------------------------------------------
% Create the grid.
g.dim = 3;
% ==========================
g.min = [0; 0; 0];
g.max = [(ysize-1)*yvoxelsize; (xsize-1)*xvoxelsize; (zsize-1)*zvoxelsize];

g.N = [ysize; xsize; zsize];

g.bdry = @addGhostExtrapolate;

g = processGrid(g);
%--------------------------------------------------------------------------
% binary circle
% center = [40;40;20];
center_binary = center_binary + g.min;
% radius = [14];
x = shapeSphere(g, center_binary, radius_binary);
index = find(x<=0);
binary_x = zeros(size(x));
binary_x(index) = 1;

SaveAnalyze(uint32(binary_x), header, 'x.hdr', 'Grey'); 
%--------------------------------------------------------------------------
% signed pressure force
SignedPressureForce = zeros(size(x));

index = find(x==0);
SignedPressureForce(index) = 0;

index = find(x>0);
SignedPressureForce(index) = -1;

index = find(x<0);
SignedPressureForce(index) = 1;

SaveAnalyze(SignedPressureForce, header, 'SignedPressureForce.hdr', 'Real'); 

%--------------------------------------------------------------------------
% create SDF
% center_sdf = [22;22;22];
center_sdf = center_sdf + g.min;
% radius_sdf = [16];
SDF = shapeSphere(g, center_sdf, radius_sdf);
data0 = SDF;

SaveAnalyze(data0, header, 'data0.hdr', 'Real'); 

%--------------------------------------------------------------------------
normalSpeed = SignedPressureForce;
bValue = 0.1;

% perform the level set propagation

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

schemeData.header = header;
schemeData.resultDir = resultDir;

schemeData = processTopology(schemeData);

%---------------------------------------------------------------------------
% Let the integrator know what function to call.
integratorOptions = odeCFLset(integratorOptions, ...
                              'postTimestep', @CheckTopology);
%--------------------------------------------------------------------------
% Initialize Display
f = figure;
axis square;  axis manual;
axis([1 g.shape(1) 1 g.shape(2) 1 g.shape(3)]);

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
    figure(f);
%     figure(f);
    
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
    view(figureView);
    drawnow;
    
    filename = ['data' num2str(tNow) '.hdr'];
    SaveAnalyze(data, header, filename, 'Real'); 
    
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