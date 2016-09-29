function [frame_index, dice, FP, FN, BSE] = ComputeMoCoPerformance_Perfusion_KLT_evaluation(dataDir, grappa, roi, xvoxelsize, yvoxelsize, withMoCo)
% Compute the moco performance statistics

%% epi

% --------------------------

slc0 = analyze75read(fullfile(dataDir, grappa, 'input_for_moco_0_MAG.hdr'));

if(withMoCo)
    dx = analyze75read(fullfile(dataDir, grappa, 'deformation_field_x_0.hdr'));
    dy = analyze75read(fullfile(dataDir, grappa, 'deformation_field_y_0.hdr'));
end

slc0_epi = analyze75read(fullfile(dataDir, roi, 'linput_for_moco_0_MAG_epi.hdr'));   
slc0_epi_segSDF = ComputeApproximatedSDF_2D(slc0_epi); 

if(withMoCo)
    slc0_epi_segSDF_moco = Matlab_gt_apply_deformation_field_reg_2D_series(single(slc0_epi_segSDF), single(dx), single(dy));    
else
    slc0_epi_segSDF_moco = slc0_epi_segSDF;
end

[frame_index_slc0_epi, dice_slc0_epi, FP_slc0_epi, FN_slc0_epi, BSE_slc0_epi] = ComputeMoCoPerformance_SDF_method_2D(slc0_epi, slc0_epi_segSDF_moco, [xvoxelsize yvoxelsize]);

% --------------------------

slc1 = analyze75read(fullfile(dataDir, grappa, 'input_for_moco_1_MAG.hdr'));
if(withMoCo)
    dx = analyze75read(fullfile(dataDir, grappa, 'deformation_field_x_1.hdr'));
    dy = analyze75read(fullfile(dataDir, grappa, 'deformation_field_y_1.hdr'));
end

slc1_epi = analyze75read(fullfile(dataDir, roi, 'linput_for_moco_1_MAG_epi.hdr'));   
slc1_epi_segSDF = ComputeApproximatedSDF_2D(slc1_epi);      

if(withMoCo)
    slc1_epi_segSDF_moco = Matlab_gt_apply_deformation_field_reg_2D_series(single(slc1_epi_segSDF), single(dx), single(dy));    
else
    slc1_epi_segSDF_moco = slc1_epi_segSDF;
end

[frame_index_slc1_epi, dice_slc1_epi, FP_slc1_epi, FN_slc1_epi, BSE_slc1_epi] = ComputeMoCoPerformance_SDF_method_2D(slc1_epi, slc1_epi_segSDF_moco, [xvoxelsize yvoxelsize]);

% --------------------------

slc2 = analyze75read(fullfile(dataDir, grappa, 'input_for_moco_2_MAG.hdr'));
if(withMoCo)
    dx = analyze75read(fullfile(dataDir, grappa, 'deformation_field_x_2.hdr'));
    dy = analyze75read(fullfile(dataDir, grappa, 'deformation_field_y_2.hdr'));
end

slc2_epi = analyze75read(fullfile(dataDir, roi, 'linput_for_moco_2_MAG_epi.hdr'));   
slc2_epi_segSDF = ComputeApproximatedSDF_2D(slc2_epi);      
if(withMoCo)
    slc2_epi_segSDF_moco = Matlab_gt_apply_deformation_field_reg_2D_series(single(slc2_epi_segSDF), single(dx), single(dy));    
else
    slc2_epi_segSDF_moco = slc2_epi_segSDF;
end

[frame_index_slc2_epi, dice_slc2_epi, FP_slc2_epi, FN_slc2_epi, BSE_slc2_epi] = ComputeMoCoPerformance_SDF_method_2D(slc2_epi, slc2_epi_segSDF_moco, [xvoxelsize yvoxelsize]);

%% endo

% ---------------------------------

slc0 = analyze75read(fullfile(dataDir, grappa, 'input_for_moco_0_MAG.hdr'));

if(withMoCo)
    dx = analyze75read(fullfile(dataDir, grappa, 'deformation_field_x_0.hdr'));
    dy = analyze75read(fullfile(dataDir, grappa, 'deformation_field_y_0.hdr'));
end

slc0_endo = analyze75read(fullfile(dataDir, roi, 'linput_for_moco_0_MAG_endo.hdr'));   
slc0_endo_segSDF = ComputeApproximatedSDF_2D(slc0_endo);      

slc0_endo_segSDF_moco = slc0_endo_segSDF;

if(withMoCo)
    slc0_endo_segSDF_moco = Matlab_gt_apply_deformation_field_reg_2D_series(single(slc0_endo_segSDF), single(dx), single(dy));    
end

[frame_index_slc0_endo, dice_slc0_endo, FP_slc0_endo, FN_slc0_endo, BSE_slc0_endo] = ComputeMoCoPerformance_SDF_method_2D(slc0_endo, slc0_endo_segSDF_moco, [xvoxelsize yvoxelsize]);

% ---------------------------------

slc1 = analyze75read(fullfile(dataDir, grappa, 'input_for_moco_1_MAG.hdr'));

if(withMoCo)
    dx = analyze75read(fullfile(dataDir, grappa, 'deformation_field_x_1.hdr'));
    dy = analyze75read(fullfile(dataDir, grappa, 'deformation_field_y_1.hdr'));
end

slc1_endo = analyze75read(fullfile(dataDir, roi, 'linput_for_moco_1_MAG_endo.hdr'));   
slc1_endo_segSDF = ComputeApproximatedSDF_2D(slc1_endo);      

slc1_endo_segSDF_moco = slc1_endo_segSDF;
if(withMoCo)
    slc1_endo_segSDF_moco = Matlab_gt_apply_deformation_field_reg_2D_series(single(slc1_endo_segSDF), single(dx), single(dy));    
end

[frame_index_slc1_endo, dice_slc1_endo, FP_slc1_endo, FN_slc1_endo, BSE_slc1_endo] = ComputeMoCoPerformance_SDF_method_2D(slc1_endo, slc1_endo_segSDF_moco, [xvoxelsize yvoxelsize]);

% ---------------------------------

slc2 = analyze75read(fullfile(dataDir, grappa, 'input_for_moco_2_MAG.hdr'));

if(withMoCo)
    dx = analyze75read(fullfile(dataDir, grappa, 'deformation_field_x_2.hdr'));
    dy = analyze75read(fullfile(dataDir, grappa, 'deformation_field_y_2.hdr'));
end

slc2_endo = analyze75read(fullfile(dataDir, roi, 'linput_for_moco_2_MAG_endo.hdr'));   
slc2_endo_segSDF = ComputeApproximatedSDF_2D(slc2_endo);      

slc2_endo_segSDF_moco = slc2_endo_segSDF;
if(withMoCo)
    slc2_endo_segSDF_moco = Matlab_gt_apply_deformation_field_reg_2D_series(single(slc2_endo_segSDF), single(dx), single(dy));    
end

[frame_index_slc2_endo, dice_slc2_endo, FP_slc2_endo, FN_slc2_endo, BSE_slc2_endo] = ComputeMoCoPerformance_SDF_method_2D(slc2_endo, slc2_endo_segSDF_moco, [xvoxelsize yvoxelsize]);

% --------------------------------

frame_index = cell(3, 2);
frame_index{1,1} = frame_index_slc0_epi;
frame_index{2,1} = frame_index_slc1_epi;
frame_index{3,1} = frame_index_slc2_epi;
frame_index{1,2} = frame_index_slc0_endo;
frame_index{2,2} = frame_index_slc1_endo;
frame_index{3,2} = frame_index_slc2_endo;

dice = cell(3, 2);
dice{1,1} = dice_slc0_epi;
dice{2,1} = dice_slc1_epi;
dice{3,1} = dice_slc2_epi;
dice{1,2} = dice_slc0_endo;
dice{2,2} = dice_slc1_endo;
dice{3,2} = dice_slc2_endo;

FP = cell(3, 2);
FP{1,1} = FP_slc0_epi;
FP{2,1} = FP_slc1_epi;
FP{3,1} = FP_slc2_epi;
FP{1,2} = FP_slc0_endo;
FP{2,2} = FP_slc1_endo;
FP{3,2} = FP_slc2_endo;

FN = cell(3, 2);
FN{1,1} = FN_slc0_epi;
FN{2,1} = FN_slc1_epi;
FN{3,1} = FN_slc2_epi;
FN{1,2} = FN_slc0_endo;
FN{2,2} = FN_slc1_endo;
FN{3,2} = FN_slc2_endo;

BSE = cell(3, 2);
BSE{1,1} = BSE_slc0_epi;
BSE{2,1} = BSE_slc1_epi;
BSE{3,1} = BSE_slc2_epi;
BSE{1,2} = BSE_slc0_endo;
BSE{2,2} = BSE_slc1_endo;
BSE{3,2} = BSE_slc2_endo;

