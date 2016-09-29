function [res_Myo, res_BL, res_PD] = ComputeMoCoPerformance_Perfusion_Slice_evaluation(slc, dataDir, grappa, roi, start_FP, end_FP, xvoxelsize, yvoxelsize, withMoCo, recomputeSDF)
% Compute the moco performance statistics
% start and end indexes for first pass

%% deformation field

data = analyze75read(fullfile(dataDir, grappa, ['input_for_moco_' num2str(slc) '_MAG.hdr']));

if(withMoCo)
    dx = analyze75read(fullfile(dataDir, grappa, ['deformation_field_x_' num2str(slc) '.hdr']));
    dy = analyze75read(fullfile(dataDir, grappa, ['deformation_field_y_' num2str(slc) '.hdr']));
end

% --------------------------
% epi

epi = fullfile(dataDir, roi, ['linput_for_moco_' num2str(slc) '_MAG_epi']);
slc_epi = analyze75read([epi '.hdr']); slc_epi = double(slc_epi);   
frame_index_epi = ComputeMoCoPerformance_find_seg_frames(slc_epi);

filename = [epi '.mat'];
if (isFileExist(filename) & ~recomputeSDF)
    load(filename);
else
    slc_epi_segSDF = ComputeApproximatedSDF_2D(slc_epi); 
    save(filename, 'slc_epi_segSDF');
end

if(withMoCo)
    slc_epi_segSDF_moco = Matlab_gt_apply_deformation_field_reg_2D_series(single(slc_epi_segSDF), single(dx), single(dy));    
else
    slc_epi_segSDF_moco = slc_epi_segSDF;
end

% --------------------------
% endo

endo = fullfile(dataDir, roi, ['linput_for_moco_' num2str(slc) '_MAG_endo']);

slc_endo = analyze75read([endo '.hdr']); slc_endo = double(slc_endo);  
frame_index_endo = ComputeMoCoPerformance_find_seg_frames(slc_endo);

filename = [endo '.mat'];
if (isFileExist(filename) & ~recomputeSDF)
    load(filename);
else
    slc_endo_segSDF = ComputeApproximatedSDF_2D(slc_endo);      
    save(filename, 'slc_endo_segSDF');
end

slc_endo_segSDF_moco = slc_endo_segSDF;

if(withMoCo)
    slc_endo_segSDF_moco = Matlab_gt_apply_deformation_field_reg_2D_series(single(slc_endo_segSDF), single(dx), single(dy));    
end

% --------------------------
% myocardium

[frame_index_common, ind_epi, ind_endo] = intersect(frame_index_epi, frame_index_endo);

slc_endo_epi = zeros(size(slc_epi));

for n=1:numel(frame_index_common)
    slc_endo_epi(:,:,frame_index_common(n)) = slc_endo(:,:,frame_index_common(n)) + slc_epi(:,:,frame_index_common(n));
end

slc_endo_epi(find(slc_endo_epi(:)>0)) = 100;

for n=1:numel(frame_index_common)
    slc_endo_epi(:,:,frame_index_common(n)) = slc_endo_epi(:,:,frame_index_common(n)) - slc_endo(:,:,frame_index_common(n));
end

myo = fullfile(dataDir, roi, ['linput_for_moco_' num2str(slc) '_MAG_myo']);

filename = [myo '.mat'];
if (isFileExist(filename) & ~recomputeSDF)
    load(filename);
else
    slc_endo_epi_segSDF = ComputeApproximatedSDF_2D(slc_endo_epi);
    save(filename, 'slc_endo_epi_segSDF');
end

if(withMoCo)
    slc_endo_epi_segSDF_moco = Matlab_gt_apply_deformation_field_reg_2D_series(single(slc_endo_epi_segSDF), single(dx), single(dy));    
else
    slc_endo_epi_segSDF_moco = slc_endo_epi_segSDF;
end

[frame_index_slc_endo_epi, dice_slc_endo_epi, FP_slc_endo_epi, FN_slc_endo_epi, BSE_slc_endo_epi] = ComputeMoCoPerformance_SDF_method_2D(slc_endo_epi, slc_endo_epi_segSDF_moco, [xvoxelsize yvoxelsize]);

% --------------------------
% Base line

BL = fullfile(dataDir, roi, ['linput_for_moco_' num2str(slc) '_MAG_BL']);

slc_BL = analyze75read([BL '.hdr']); 

filename = [BL '.mat'];
if (isFileExist(filename) & ~recomputeSDF)
    load(filename);
else
    slc_BL_segSDF = ComputeApproximatedSDF_2D(slc_BL);
    save(filename, 'slc_BL_segSDF');
end

if(withMoCo)
    slc_BL_segSDF_moco = Matlab_gt_apply_deformation_field_reg_2D_series(single(slc_BL_segSDF), single(dx), single(dy));    
else
    slc_BL_segSDF_moco = slc_BL_segSDF;
end

[frame_index_slc_BL, dice_slc_BL, FP_slc_BL, FN_slc_BL, BSE_slc_BL] = ComputeMoCoPerformance_SDF_method_2D(slc_BL, slc_BL_segSDF_moco, [xvoxelsize yvoxelsize]);

%% PD

PD = fullfile(dataDir, roi, ['linput_for_moco_' num2str(slc) '_MAG_PD']);

slc_PD = analyze75read([PD '.hdr']);   

filename = [PD '.mat'];
if (isFileExist(filename) & ~recomputeSDF)
    load(filename);
else
    slc_PD_segSDF = ComputeApproximatedSDF_2D(slc_PD);
    save(filename, 'slc_PD_segSDF');
end

if(withMoCo)
    slc_PD_segSDF_moco = Matlab_gt_apply_deformation_field_reg_2D_series(single(slc_PD_segSDF), single(dx), single(dy));    
else
    slc_PD_segSDF_moco = slc_PD_segSDF;
end

[frame_index_slc_PD, dice_slc_PD, FP_slc_PD, FN_slc_PD, BSE_slc_PD] = ComputeMoCoPerformance_SDF_method_2D(slc_PD, slc_PD_segSDF_moco, [xvoxelsize yvoxelsize]);

%% get results
res_Myo = cell(5, 1);
res_Myo{1} = frame_index_slc_endo_epi;
res_Myo{2} = dice_slc_endo_epi;
res_Myo{3} = FP_slc_endo_epi;
res_Myo{4} = FN_slc_endo_epi;
res_Myo{5} = BSE_slc_endo_epi;

res_BL = cell(5, 1);
res_BL{1} = frame_index_slc_BL;
res_BL{2} = dice_slc_BL;
res_BL{3} = FP_slc_BL;
res_BL{4} = FN_slc_BL;
res_BL{5} = BSE_slc_BL;

res_PD = cell(5, 1);
res_PD{1} = frame_index_slc_PD;
res_PD{2} = dice_slc_PD;
res_PD{3} = FP_slc_PD;
res_PD{4} = FN_slc_PD;
res_PD{5} = BSE_slc_PD;


