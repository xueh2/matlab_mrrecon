
function option = CreateDynamicReconOption()
% create the default option for the dynamic recon
% options : the options to control the behaviors
% options.reconStrategy = 'TGRAPPA' or 'TSENSE' or 'HTGRAPPA' or 'IcePAT'
% options.sensitivity = [] or pre-scan kspace for sensitivity [Nfe Npe numOfCoil]
% options.reference = 'AverageAll' or 'SlidingWindow'
% options.subtractMean = 1 if the subtract mean recon is performed
% options.timing = 1 if the recon time is measured
% options.preWhitening = 1 if performing the noise prewhitening
% options.zeroFilledSize = [newNfe newNpe], non-empty if the image will be resized using zero filling
%
% For 'TGRAPPA' or 'HTGRAPPA'
% ----------------------------------------
% options.linesUsedForGRAPPA = 'All' or 'frame'
% Lines in the reference images used to compute the grappa kernels, 'All' means all data lines are used; 'frame' means for a 
% frame, only lines that were really collected are used
% options.updateRef = 'Once', or 'Most' 'Once' means to compute the grappa kernel only once from reference image
% 'Most' means to compute the grappa kernel for every frame
%
% options.minFEUsed, options.maxFEUsed : the FE points used to estimated grappa kernels are between [options.minFEUsed options.maxFEUsed]
% options.blockCoord : line blocks used to estimated grappa kernels
% options.missingLines : missed lines in a block
% options.firstDataLineInBlcok : index of first data line in the block
% options.kernalIndFE : the offset of 2D grappa kernels along the FE direction
% options.useACSLines : 1 if the acs lines used to replace the estimated kspace lines
% options.regularization : 'svd' or 'tikhonov'
% options.thresReg : threshold for regulization
% options.overDetermineRatio : the overdetermine ratio = number of equantions / number of coefficients for grappa 
% options.kspaceCenterFE : the kspace centre along the FE direction
% options.paritalFourierFE : 1 if partial fourier is used along the FE direction
% options.phaseEncodingUsed : the portion of PE lines used, [0 1], 1 means 100% lines used, 0.5 means 50% lines used
%
% options.rawFilterFE : 'None, Unity, Gaussian, Hanning, Hamming, Tukey, Sinus'
% options.rawFilterFEStrength : 'Weak, Medium, Strong, VeryStrong'
% options.rawFilterPE : 'None, Unity, Gaussian, Hanning, Hamming, Tukey, Sinus'
% options.rawFilterPEStrength : 'Weak, Medium, Strong, VeryStrong'
%
% options.ascFilterFE : 'None, Unity, Gaussian, Hanning, Hamming, Tukey, Sinus'
% options.ascFilterFEStrength : 'Weak, Medium, Strong, VeryStrong'
% options.ascFilterPE : 'None, Unity, Gaussian, Hanning, Hamming, Tukey, Sinus'
% options.ascFilterPEStrength : 'Weak, Medium, Strong, VeryStrong'
%
% For 'SlidingWindow'
% ----------------------------------------
% options.numOfBlocksForRef = number of frame blocks used to estimated reference images if options.reference == 'SlidingWindow'
% A frame block is defined as a set of consecutive reductionFactor frames
% For 'SlidingWindow' & 'TSENSE', the sensitivity map is updated whenever a new block of frames is added or dropped
% For 'SlidingWindow' & 'TGRAPPA' or 'HTGRAPPA', the grappa kernel is updated whenever a new block of frames is added or dropped

option = struct('reconStrategy', 'TGRAPPA', 'sensitivity', [], 'reference', 'AverageAll', 'subtractMean', 1, 'preWhitening', 0, ...
    'linesUsedForGRAPPA', 'frame', 'updateRef', 'Once', 'minFEUsed', 0.1, 'maxFEUsed', 0.9, 'blockCoord', [], 'missingLines', [], ...
    'firstDataLineInBlcok', [], 'kernalIndFE', [2:-1:-2], 'useACSLines', 1, 'regularization', 'svd', 'thresReg', 1e-3, ...
    'overDetermineRatio', 10, 'kspaceCenterFE', -1, 'paritalFourierFE', 0, 'phaseEncodingUsed', 1.0, 'zeroFilledSize', [], ...
    'rawFilterFE', 'Sinus', 'rawFilterFEStrength', 'weak', 'rawFilterPE', 'Sinus', 'rawFilterPEStrength', 'weak', ... 
    'ascFilterFE', 'None', 'ascFilterFEStrength', 'weak', 'ascFilterPE', 'None', 'ascFilterPEStrength', 'None', ... 
    'numOfBlocksForRef', 3, 'timing', 0, 'voxelsize', [1 1 1], 'iters', [32 16 8], 'sigma', 128, ...  
    'maxIterSPIRiT', 70, 'objTollSPIRiT', 0.01, 'wavWeightSPIRiT', 0.0015, 'cgSPIRiTlambda', 1e-4, 'TVWeightSPIRiT', 1e-4, 'dataWeightSPIRiT', 0.05, ... 
    'showIterSPIRiT', 0, 'stopThresSPIRiT', 1e-4, 'pocsFlagSPIRiT', 0, 'centre', 250, 'width', 500, ... 
    'continuationStep', 10, 'wavThresRatio', 2, 'KLTSen', 1, 'numOfModesKeptKLTSen', 3, 'unwarpMethod', 'fft', ... 
    'kspaceFullRef', [], 'predictionKSpace', []);
