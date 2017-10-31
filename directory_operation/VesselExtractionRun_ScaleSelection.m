
function VesselExtractionRun_ScaleSelection(currentDir, tofResampled, extractionResults)
% perform the vessel extraction

path_abo = fullfile(currentDir, tofResampled, '*.hdr' );
indir = dir(path_abo);

num_indir = length(indir);

if ( num_indir == 0 )
    disp('resampled hdr file is not available...');
    return
end

% load image

filename = fullfile(currentDir, tofResampled, indir(1).name);
[data, header] = LoadAnalyze(filename, 'Grey');

[xy, yz, zx]=getMIP3(data, header);
areaPercentage = 0.7;
background = 40;
[thresholdxy, segmentedMIPxy] = OptimalThreshold(xy, areaPercentage, background);
[thresholdyz, segmentedMIPyz] = OptimalThreshold(yz, areaPercentage, background);
[thresholdzx, segmentedMIPzx] = OptimalThreshold(zx, areaPercentage, background);
t = min([thresholdxy, thresholdyz, thresholdzx]);
t = t * 0.8;

place = find(currentDir(length(currentDir):-1:1) == filesep );
prefix = currentDir(length(currentDir)-place(1)+2:length(currentDir));

% extracting vessels
data = single(data);
maxData = max(max(max(abs(data))));
data = normalizeImage(data);

% global definition and declaration
global vasculature
global Aylward_coefficient
% vasculature
vasculature = struct('count',0, 'vessels', [], 'totalnumofpoints', 0);
onevessel = struct('ridgepoints',[], 'radius',[], 'samplepositions', [], 'len', 0, 'dist', [],...
    'points_eigenvalue', [], 'points_eigenvector1', [], 'points_eigenvector2', [], ...
    'points_eigenvector3', [], 'points_sigma', []);
vasculatureNodes = struct('count',0, 'nodes', [], 'totalnumofpoints', 0);

suffix = prefix;

seedfactor = 0.08;
samplestepnew = 0.05;

samplestep = 0.2;
numofkernels = 5;
%-------------------------------------------------------------------------%
parameters = struct('tolerance',1.0e-4, 'symmetry', 1, 'minSigma', 0.2, 'maxSigma', 3.0, 'sigma0', 0.2, 'sigmastep', 0.15,...
    'beta', 0.2, 'N', 180, 'stepsizethreshold', 0.001, 'measuredstep', 5, 'tk_straightline_low', 0.7,...
    'tk_straightline_high', 0.9, 'unusededgewidth', 2, 'maxStep', 1, 'recomputescale', 10,...
    'minimalpercent', 1, 'scalepercent', 0.8, 'function', @ridgetranverse_direct1voxel9_ssdofHM_adaptive_nolinearsearch);
                    
radiusparameters = struct('minSigma', 0.15, 'maxSigma', 3.0, 'sigmastep', 0.05, ...
    'N', 360, 'smoothfactor', 0.8, 'minimalpercent', 0.75, 'samplestep', samplestep, 'numofkernels', numofkernels);
factor = 2;

doradius = 0;

% presetScales = [0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2 2.0 2.2 2.4 2.6];

presetScales = [ 0.4 0.8 1.2 1.6 2.0 2.4];

t = t / maxData;

background = t * 0.8;

zbsflag = 1;
zbs_usedSeedsflag = 1;
coefficientflag = 1;
coefficientblurredflag = 1;

zbs_leftup = [];
zbs_rightdown = [];

lowest_seeds = 1000;
highest_seeds = 3000;

[vasculature, coefficient, vasculatureDetected] = RidgeDetection5_ScaleSelection(data, header, zbs_leftup, zbs_rightdown,...
    seedfactor, samplestepnew, parameters, radiusparameters, factor, suffix,...
    doradius, zbsflag, zbs_usedSeedsflag, coefficientflag, coefficientblurredflag, background, t, presetScales, lowest_seeds, highest_seeds);

% save vasculature into the extractionResults directory...
vesselfile = fullfile(currentDir, extractionResults, 'vasculature_scaleselection2.mat');
save(vesselfile, 'vasculature');
return;