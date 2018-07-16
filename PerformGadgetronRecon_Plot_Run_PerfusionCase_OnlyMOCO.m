
function PerformGadgetronRecon_Plot_Run_PerfusionCase_OnlyMOCO(perf_cases, rest_cases, resDir, figDir, only_reviewing)
% PerformGadgetronRecon_Plot_Run_PerfusionCase_OnlyMOCO(perf_cases, rest_cases, resDir, figDir, only_reviewing)

if(nargin < 5)
    only_reviewing = 0;
end

saveDir = fullfile('\\hl-share\RawMRI\Lab-Kellman\Share', figDir);
if(isunix())
    saveDir = fullfile('/mnt/Lab-Kellman/Share', figDir);
end
mkdir(saveDir)

cases = [];

startN = 1
endN = size(perf_cases, 1)
for ii=startN:endN
    cases = [cases; perf_cases(ii,2); perf_cases(ii,3)];
end

endN = size(rest_cases, 1)
for ii=startN:endN
    cases = [cases; rest_cases(ii,2)];
end

startN = 1
endN = size(cases, 1)
for ii=startN:endN
    disp([num2str(ii-startN+1) ' out of ' num2str(endN-startN+1) ' - ' cases(ii)]);

    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(cases{ii});
    
    debug_dir = fullfile(resDir, study_dates, cases{ii}, 'DebugOutput');
    cd(debug_dir)
    s1 = analyze75read('input_for_moco_0_MAG');
    s2 = analyze75read('input_for_moco_1_MAG');
    s3 = analyze75read('input_for_moco_2_MAG');
    
    f1 = analyze75read('flow_map_0_MAG'); f1 = f1(:,:,end);
    f2 = analyze75read('flow_map_1_MAG'); f2 = f2(:,:,end);
    f3 = analyze75read('flow_map_2_MAG'); f3 = f3(:,:,end);

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

    %% try old moco
    num_levels_moco_among_model = 5; 
    regularization_scale_factor_moco_among_model = 1.5; 
    dissimilarity_LocalCCR_sigmaArg_moco_among_model = 2; 

    tic
    [dx1_old, dy1_old, moco1_old, dxInv1_old, dyInv1_old] = Matlab_gt_perfusion_model_moco(single(s1), keyframe, level, max_iter_num_pyramid_level, ... 
        LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, DivergenceFreeReg, ... 
        KLT_regularization_hilbert_strength, KLT_regularization_minimal_hilbert_strength, num_levels_moco_among_model, regularization_scale_factor_moco_among_model, ... 
        dissimilarity_LocalCCR_sigmaArg_moco_among_model, numPD, DebugFolder, verbose);
    toc

    tic
    [dx2_old, dy2_old, moco2_old, dxInv2_old, dyInv2_old] = Matlab_gt_perfusion_model_moco(single(s2), keyframe, level, max_iter_num_pyramid_level, ... 
        LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, DivergenceFreeReg, ... 
        KLT_regularization_hilbert_strength, KLT_regularization_minimal_hilbert_strength, num_levels_moco_among_model, regularization_scale_factor_moco_among_model, ... 
        dissimilarity_LocalCCR_sigmaArg_moco_among_model, numPD, DebugFolder, verbose);
    toc

    tic
    [dx3_old, dy3_old, moco3_old, dxInv3_old, dyInv3_old] = Matlab_gt_perfusion_model_moco(single(s3), keyframe, level, max_iter_num_pyramid_level, ... 
        LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, DivergenceFreeReg, ... 
        KLT_regularization_hilbert_strength, KLT_regularization_minimal_hilbert_strength, num_levels_moco_among_model, regularization_scale_factor_moco_among_model, ... 
        dissimilarity_LocalCCR_sigmaArg_moco_among_model, numPD, DebugFolder, verbose);
    toc
    
    %% make figures
    h = figure('Name',cases{ii},'NumberTitle','off'); 
    res = cat(4, s1(:,:,numPD+1:end), s2(:,:,numPD+1:end), s3(:,:,numPD+1:end), moco1(:,:,numPD+1:end), moco2(:,:,numPD+1:end), moco3(:,:,numPD+1:end), moco1_old(:,:,numPD+1:end), moco2_old(:,:,numPD+1:end), moco3_old(:,:,numPD+1:end));
    res = flipdim(res, 1);
    res = flipdim(res, 2);
    imagescn(res, [0 250], [3 3], [8], 3);
    
    fmap = cat(3, f1, f2, f3);
    fmap = flipdim(fmap, 1);
    fmap = flipdim(fmap, 2);
    h2 = figure('Name', [cases{ii} '_FlowMap'],'NumberTitle','off'); 
    imagescn(fmap, [0 8], [1 3], [8]); PerfColorMap;
    
    savefig(h, fullfile(saveDir, cases{ii}));        
    savefig(h2, fullfile(saveDir, [cases{ii} '_FlowMap'] ));
    
    tic
    save perf_moco s1 s2 s3 moco1 moco2 moco3 dx1 dx2 dx3 dy1 dy2 dy3 dxInv1 dxInv2 dxInv3 dyInv1 dyInv2 dyInv3
    toc
    
    if(only_reviewing)
        pause;
    end
    closeall
    
    disp(['====================================================================================================']);
end
