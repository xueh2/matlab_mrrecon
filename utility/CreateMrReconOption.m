function option = CreateMrReconOption(FOV_reduction_factor, FOV_reduction_factor_PAR, coilMapMethod, dstChaThres, performSENSE, performGRAPPA, performVICPAAS, SamplesInScan, KSpaceCentreColumn, MaxKSpaceLineNo, KSpaceCentreLineNo, KSpaceCentrePartitionNo)
% create the default option structure for MrRecon

option.KernelSize = [5 4 4]
option.thresReg  = 1e-4
option.KernelPattern = [-FOV_reduction_factor:FOV_reduction_factor:2*FOV_reduction_factor];
option.KernelPatternPAR = [-FOV_reduction_factor_PAR:FOV_reduction_factor_PAR:2*FOV_reduction_factor_PAR];
option.OutPattern = [0:FOV_reduction_factor-1];
option.OutPatternPAR = [0:FOV_reduction_factor_PAR-1];
option.GrappaOnly = 0;
option.useAveAllKernel = 1;

% v-spirit part
option.KernelSizeSPIRiT = [5 5 5];
option.OutKernelSPIRiT = [3 3 3];

option.KernelSizeSPIRiT_TPAT = [7 7 5];
option.OutKernelSPIRiT_TPAT = [1 1 1];

option.maxIterSPIRiT = 150;
option.stopThresSPIRiT = 1e-5;
option.pocsFlagSPIRiT = 0;
option.continuationStep = 10;
option.wavThresRatio = 1;
option.unwarpMethod = 'fft_AllGPU';

option.thresRegSPIRiT = 0.005;
option.thresRegSPIRiT_TPAT = 0.005;

option.wavWeightSPIRiT = 0.001
option.TVWeightSPIRiT = 0
option.dataWeightSPIRiT = -1
option.temporalScalingFactorSPIRiT = 5; 

option.Itnlim = 10;
option.objTollSPIRiT = 0.1;

option.multiTaskRecon = 1;
option.performReconWithMoCo = 0;
option.UseCoilSensitivity = 1;
option.computeSPIRITLinear = 1;
option.performCenteredFFT = 1;
option.TwoDPlusTRecon = 1;
option.printInfo = 0;
option.tspirit = 0;
option.dynamicKernel = 1;
option.vicpaa = 1;
option.lsqrForLinear = 0;

% common parameters
option.csmMethod = coilMapMethod;
option.dstCha  = -1;
option.dstChaThres  = dstChaThres;

option.rawFilterFE = 'Gaussian';
option.rawFilterFEStrength = 'Weak';
option.rawFilterPE = 'Gaussian';
option.rawFilterPEStrength = 'Weak';
option.rawFilterPAR = 'Gaussian';
option.rawFilterPARStrength = 'Weak';

option.acsFilterFE = 'Hanning';
option.acsFilterFEStrength = 'Strong';
option.acsFilterPE = 'Hanning';
option.acsFilterPEStrength = 'Strong';
option.acsFilterPAR = 'Hanning';
option.acsFilterPARStrength = 'Strong';

option.KLTSen = 1;
option.numOfModesKeptKLTSen = 3;
option.spatialSmoothingKernel  = 5;
option.kspaceCenterFE = KSpaceCentreColumn;

option.kSizeEigenVectorCoilSensitivity = [5 5];
option.percentageEigenVectorCoilSensitivity = 96;
option.thresForSecondSenMap = 0.9;

option.zeroFilledSize = []

option.readoutDecoupled = 1;

option.centre = 1024;
option.width = 1024;    

option.performSENSE = performSENSE;
option.performGRAPPA = performGRAPPA;
option.performVICPAAS = performVICPAAS;

option.SamplesInScan = SamplesInScan;
option.KSpaceCentreColumn = KSpaceCentreColumn;
option.MaxKSpaceLineNo = MaxKSpaceLineNo;
option.KSpaceCentreLineNo = KSpaceCentreLineNo;
option.KSpaceCentrePartitionNo = KSpaceCentrePartitionNo;

option.handleAsymmetricEcho = 1;
option.handlePartialFourierPE = 1;
option.handlePartialFourierPAR = 1;

% for the sense recon
option.niter = 30;
option.lambda = 1e-6;
option.alpha = 0.1;
option.senseSubMean = 1;
option.thres = 1e-2/2;

% for the partial fourier handling
% option.PFMethod = 'FengHuang';
% option.PFMethod = 'SPIRiTPFCorr';
option.PFMethod = 'None';

option.homodyne.iters = 6;
option.homodyne.thres = 1e-3;

option.pocs.iters = 6;
option.pocs.thres = 1e-3;

option.FengHuang.kSize = [5 5 5];
option.FengHuang.thresReg = 0.01;

option.SPIRiTPFCorr.kSize = [5 5];
option.SPIRiTPFCorr.OutPattern = [1 1];
option.SPIRiTPFCorr.thresReg = 0.01;
option.SPIRiTPFCorr.performCenteredFFT = 1;
option.SPIRiTPFCorr.maxIter = 70;
option.SPIRiTPFCorr.stopThres = 1e-3;
option.SPIRiTPFCorr.printInfo = 0;

option
