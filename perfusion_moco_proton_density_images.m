function [moco_perf, moco_pd_klt, moco_pd_new, moco_pd_new_no_lastpd_firstperf, moco_pd_new_after_deform_check, moco_klt, moco_co_reg, mean_deform] = perfusion_moco_proton_density_images(perf, pd, roi, max_mean_deform, dx, dy)
% [moco_perf, moco_pd_klt, moco_pd_new, moco_klt, moco_co_reg] = perfusion_moco_proton_density_images(perf, pd)
% perform PD moco using the deformation field concatenation method

if(nargin<5)
    dx = [];
end

if(nargin<6)
    dy = [];
end

RO = size(perf, 1);
E1 = size(perf, 2);

ps_x = roi(1);
pe_x = roi(2);
ps_y = roi(3);
pe_y = roi(4);

keyframe = 32; 
level = 4; 
max_iter_num_pyramid_level = [100 100 100 100]; 
LocalCCR_sigmaArg = 2.0; 
BidirectionalReg = 1; 
dissimilarity_thres = 1e-6; 
DivergenceFreeReg = 0; 
KLT_regularization_hilbert_strength = [18 16 14 12 10]; 
KLT_regularization_minimal_hilbert_strength = 9; 
num_levels_moco_among_model = 5; 
regularization_scale_factor_moco_among_model = 1.7; 
dissimilarity_LocalCCR_sigmaArg_moco_among_model = 1.5; 	
DebugFolder = []; 
verbose = 0;

num_perf_for_PD = size(pd, 3); 
perf_pd = cat(3, pd, perf);

if(~isempty(dx) & ~isempty(dy))
    perf_pd_dx = zeros(size(perf_pd));
    perf_pd_dy = zeros(size(perf_pd));
    moco_perf_pd = perf_pd;
    
    a = Matlab_gt_apply_deformation_field_reg_2D_series(double(perf), dx, dy);
    
    perf_pd_dx(:,:,num_perf_for_PD+1:end) = dx;
    perf_pd_dy(:,:,num_perf_for_PD+1:end) = dy;
    
    moco_perf_pd(:,:,num_perf_for_PD+1:end) = a;
else
    tic
    [perf_pd_dx, perf_pd_dy, moco_perf_pd, perf_pd_dxInv, perf_pd_dyInv] = Matlab_gt_perfusion_model_moco(single(perf_pd), keyframe, level, max_iter_num_pyramid_level, LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, DivergenceFreeReg, KLT_regularization_hilbert_strength, KLT_regularization_minimal_hilbert_strength, num_levels_moco_among_model, regularization_scale_factor_moco_among_model, dissimilarity_LocalCCR_sigmaArg_moco_among_model, num_perf_for_PD, DebugFolder, verbose); 
    toc
end

moco_pd_klt = moco_perf_pd(:,:,1:num_perf_for_PD+1);
moco_perf = moco_perf_pd(:,:,num_perf_for_PD+1:end);

% test pd to all perf
min_deform = zeros(size(pd,3),1);
best_perf_ind = zeros(size(pd,3),1);
for tt=1:size(pd, 3)
    im = cat(3, pd(:,:,tt), perf);
    im = im(ps_x:pe_x-1,ps_y:pe_y-1,:);

    keyFrame = 0; 
    strategy = 'FixedReference'; 
    dissimilarity = 'LocalCCR'; 
    level = 4; 
    max_iter_num_pyramid_level = [32 64 100 100]; 
    regularization_hilbert_strength = 16.0; 
    LocalCCR_sigmaArg = 4.0; 
    BidirectionalReg = 1; 
    dissimilarity_thres = 1e-6; 
    div_num = 3; 
    inverse_deform_enforce_iter = 10; 
    inverse_deform_enforce_weight = 0.5; 
    DivergenceFreeReg = 0; 
    DebugFolder = []; 
    verbose = 0; 
    [dx, dy, moco_lastPd_firstPerf, dxInv, dyInv] = Matlab_gt_deformation_field_reg_2D_series(double(im), keyFrame, strategy, dissimilarity, level, max_iter_num_pyramid_level, regularization_hilbert_strength, LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, div_num, inverse_deform_enforce_iter, inverse_deform_enforce_weight, DivergenceFreeReg, DebugFolder, verbose); 

    mean_deform = zeros(size(dx, 3)-1, 1);
    for ii=2:size(dx, 3)    
        aa = sqrt(dx(:,:,ii).*dx(:,:,ii) + dy(:,:,ii).*dy(:,:,ii));;
        mean_deform(ii-1,1) = mean(aa(:));    
    end

    [min_deform(tt), best_perf_ind(tt)] = min(mean_deform(:));
end

[min_deform_best_pd, which_best_pd] = min(min_deform);

im = cat(3, pd(:,:,which_best_pd), perf(:,:,best_perf_ind(which_best_pd)));
keyFrame = 1; 
regularization_hilbert_strength = 32.0; 
LocalCCR_sigmaArg = 2.0;
DivergenceFreeReg = 1; 
[dx, dy, moco_lastPd_firstPerf2, dxInv, dyInv] = Matlab_gt_deformation_field_reg_2D_series(double(im), keyFrame, strategy, dissimilarity, level, max_iter_num_pyramid_level, regularization_hilbert_strength, LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, div_num, inverse_deform_enforce_iter, inverse_deform_enforce_weight, DivergenceFreeReg, DebugFolder, verbose); 

% first pd to best pd
if(which_best_pd~=1)
    im_pd = cat(3, pd(:,:,1), pd(:,:,which_best_pd));
    regularization_hilbert_strength = 16.0;
    LocalCCR_sigmaArg = 2.0; 
    DivergenceFreeReg = 0; 
    [dx_pd1st_pdLast, dy_pd1st_pdLast, moco__pd1st_pdLast, dxInv_pd1st_pdLast, dyInv_pd1st_pdLast] = Matlab_gt_deformation_field_reg_2D_series(double(im_pd), keyFrame, strategy, dissimilarity, level, max_iter_num_pyramid_level, regularization_hilbert_strength, LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, div_num, inverse_deform_enforce_iter, inverse_deform_enforce_weight, DivergenceFreeReg, DebugFolder, verbose); 
else
    dx_pd1st_pdLast = zeros(RO, E1, 2);
    dy_pd1st_pdLast = zeros(RO, E1, 2);
end

% conc deformation fields

dx_c = zeros(RO, E1, 4);
dy_c = zeros(RO, E1, 4);

dx_c(:,:,1:3) = cat(3, dx_pd1st_pdLast(:,:,1), dx(:,:,1), perf_pd_dx(:,:,num_perf_for_PD+best_perf_ind(which_best_pd)));
dy_c(:,:,1:3) = cat(3, dy_pd1st_pdLast(:,:,1), dy(:,:,1), perf_pd_dy(:,:,num_perf_for_PD+best_perf_ind(which_best_pd)));

[dx_out, dy_out] = Matlab_gt_concatenate_deform_fields(double(dx_c), double(dy_c), 3); 

moco_pd_new = Matlab_gt_apply_deformation_field_reg_2D_series(double(pd(:,:,1)), dx_out(:,:,1), dy_out(:,:,1));

% conc deformation field with mean deformation check
ps_x = roi(1);
pe_x = roi(2);
ps_y = roi(3);
pe_y = roi(4);
dx_out_cropped = dx_out(ps_x:pe_x-1,ps_y:pe_y-1);
dy_out_cropped = dy_out(ps_x:pe_x-1,ps_y:pe_y-1);

mag_deform = sqrt(dx_out_cropped.*dx_out_cropped + dy_out_cropped.*dy_out_cropped);
mean_deform = mean(mag_deform(:));

disp(['Mean deform within heart ROI is ' num2str(mean_deform)]);

% do not apply last_pd to first perf deformation field
dx(:,:,1) = 0;
dy(:,:,1) = 0;

if(mean_deform>=max_mean_deform)
    dx_c(:,:,1:3) = cat(3, dx_pd1st_pdLast(:,:,1), dx(:,:,1), perf_pd_dx(:,:,num_perf_for_PD+1));
    dy_c(:,:,1:3) = cat(3, dy_pd1st_pdLast(:,:,1), dy(:,:,1), perf_pd_dy(:,:,num_perf_for_PD+1));
else
    dx_c(:,:,1:3) = cat(3, dx_pd1st_pdLast(:,:,1), dx(:,:,1), perf_pd_dx(:,:,num_perf_for_PD+best_perf_ind(which_best_pd)));
    dy_c(:,:,1:3) = cat(3, dy_pd1st_pdLast(:,:,1), dy(:,:,1), perf_pd_dy(:,:,num_perf_for_PD+best_perf_ind(which_best_pd)));
end

[dx_out, dy_out] = Matlab_gt_concatenate_deform_fields(double(dx_c), double(dy_c), 3); 

moco_pd_new_no_lastpd_firstperf = Matlab_gt_apply_deformation_field_reg_2D_series(double(pd(:,:,1)), dx_out(:,:,1), dy_out(:,:,1));

moco_pd_new_after_deform_check = moco_pd_new;
if(mean_deform>=max_mean_deform)
    moco_pd_new_after_deform_check = moco_pd_new_no_lastpd_firstperf;
end

moco_klt = cat(3, moco_pd_klt(:,:,1), moco_perf);
moco_co_reg = cat(3, moco_pd_new_after_deform_check(:,:,1), moco_perf);
