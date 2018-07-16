
function [ori, moco, dx, dy] = PerformGadgetronRecon_Plot_Run_PerfusionCase_OnlyMOCO_OneCase(debugDir, only_reviewing)
% PerformGadgetronRecon_Plot_Run_PerfusionCase_OnlyMOCO_OneCase(debugDir, only_reviewing)

if(nargin < 2)
    only_reviewing = 0;
end

    
cd(debugDir)
s1 = analyze75read('input_for_moco_0_MAG');
s2 = analyze75read('input_for_moco_1_MAG');
s3 = analyze75read('input_for_moco_2_MAG');

keyframe = 20; 
level = 4; 
max_iter_num_pyramid_level = [32 100 100 100]; 
LocalCCR_sigmaArg = 2.0; 
BidirectionalReg = 1; 
dissimilarity_thres = 1e-6; 
DivergenceFreeReg = 0; 
KLT_regularization_hilbert_strength = [18 16 14 12 10]; 
KLT_regularization_minimal_hilbert_strength = 9; 
num_levels_moco_among_model = 3; 
regularization_scale_factor_moco_among_model = 1.25; 
dissimilarity_LocalCCR_sigmaArg_moco_among_model = 1.5; 
numPD = 3;
% DebugFolder = 'D:/gtuser/mrprogs/install/DebugFolder/'; 
DebugFolder = []; 
verbose = 1; 

tic
[dx1, dy1, moco1, dxInv1, dyInv1] = Matlab_gt_perfusion_model_moco(single(s1), keyframe, level, max_iter_num_pyramid_level, ... 
    LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, DivergenceFreeReg, ... 
    KLT_regularization_hilbert_strength, KLT_regularization_minimal_hilbert_strength, num_levels_moco_among_model, regularization_scale_factor_moco_among_model, ... 
    dissimilarity_LocalCCR_sigmaArg_moco_among_model, numPD, DebugFolder, verbose);
toc

tic
[dx2, dy2, moco2, dxInv2, dyInv2] = Matlab_gt_perfusion_model_moco(single(s2), keyframe, level, max_iter_num_pyramid_level, ... 
    LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, DivergenceFreeReg, ... 
    KLT_regularization_hilbert_strength, KLT_regularization_minimal_hilbert_strength, num_levels_moco_among_model, regularization_scale_factor_moco_among_model, ... 
    dissimilarity_LocalCCR_sigmaArg_moco_among_model, numPD, DebugFolder, verbose);
toc

tic
[dx3, dy3, moco3, dxInv3, dyInv3] = Matlab_gt_perfusion_model_moco(single(s3), keyframe, level, max_iter_num_pyramid_level, ... 
    LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, DivergenceFreeReg, ... 
    KLT_regularization_hilbert_strength, KLT_regularization_minimal_hilbert_strength, num_levels_moco_among_model, regularization_scale_factor_moco_among_model, ... 
    dissimilarity_LocalCCR_sigmaArg_moco_among_model, numPD, DebugFolder, verbose);
toc

ori = cat(4, s1, s2, s3);
moco = cat(4, moco1, moco2, moco3);
dx = cat(4, dx1, dx2, dx3);
dy = cat(4, dy1, dy2, dy3);

figure; imagescn(cat(4, s1(:,:,numPD+1:end), s2(:,:,numPD+1:end), s3(:,:,numPD+1:end), moco1(:,:,numPD+1:end), moco2(:,:,numPD+1:end), moco3(:,:,numPD+1:end)), [0 250], [2 4], [6], 3);

if(only_reviewing)
    pause;
end
closeall
