function [rad_strain, circ_strain, res_data] = strainWrapper(data, cine, mask, strain_Fn, phase, is_tag)

if is_tag
    dissimilarity = 'LocalCCR';
    level = 4;
    max_iter_num_pyramid_level = [200 200 200 200];
    regularization_hilbert_strength = [8 8 8 8];
    LocalCCR_sigmaArg = [2 2 2 2];
    BidirectionalReg = 1;
    dissimilarity_thres = 0;
    div_num = 3;
    inverse_deform_enforce_iter = 10;
    inverse_deform_enforce_weight = 0.5;
    DivergenceFreeReg = 0;
    DebugFolder = [];
    verbose = 0;
    RO = size(data,1);
    E1 = size(data,2);
    PHS = size(data,3);

    key_frame = phase;

    data1 = data;
    data2 = data;
    data2(:,:,1:key_frame-1) = data(:,:,2:key_frame);
    data2(:,:,key_frame+1:PHS) = data(:,:,key_frame:PHS-1);
    size(data1)
    size(data2)

    % data1 = repmat(a(:,:,10), [1 1 30]);
    % data2 = a;
    tic
    [dx, dy, warped, dxInv, dyInv] = Matlab_gt_deformation_field_reg_2D_series_pairwise(single(data2), single(data1), dissimilarity, level, max_iter_num_pyramid_level, regularization_hilbert_strength, LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, div_num, inverse_deform_enforce_iter, inverse_deform_enforce_weight, DivergenceFreeReg, DebugFolder, verbose);
    toc

    figure; imagescn(cat(4, data1, data2), [], [], [12], 4);
    figure; imagescn(cat(4, data2, warped), [], [], [12], 4);

    tic
    [dx_out, dy_out] = Matlab_gt_concatenate_deform_fields(double(dx), double(dy), key_frame-1);
    toc

    dx = dx_out;
    dy = dy_out;

    res_data = Matlab_gt_apply_deformation_field_reg_2D_series(double(data), dx_out, dy_out);
    res_mask = Matlab_gt_apply_deformation_field_reg_2D_series(double(mask), dx_out, dy_out);
    res_cine = Matlab_gt_apply_deformation_field_reg_2D_series(double(cine), dx_out, dy_out);
    
    [rad_strain, circ_strain] = strain_Fn(mask(:, :, phase), dx, dy);
else
    phase = 10;
    dissimilarity = 'LocalCCR';
    level = 4;
    max_iter_num_pyramid_level = [200 200 200 200];
    regularization_hilbert_strength = [16 16 16 16];
    LocalCCR_sigmaArg = [3 3 3 3];
    BidirectionalReg = 1;
    dissimilarity_thres = 0;
    div_num = 3;
    inverse_deform_enforce_iter = 10;
    inverse_deform_enforce_weight = 0.5;
    DivergenceFreeReg = 0;
    DebugFolder = [];
    verbose = 0;
    RO = size(data,1);
    E1 = size(data,2);
    PHS = size(data,3);

    data1 = repmat(data(:,:,phase), [1 1 PHS]);
    data2 = data;
    
    size(data1)
    size(data2)
    tic
    [dx, dy, warped, ~, ~] = Matlab_gt_deformation_field_reg_2D_series_pairwise(single(data1), single(data2), dissimilarity, level, max_iter_num_pyramid_level, regularization_hilbert_strength, LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, div_num, inverse_deform_enforce_iter, inverse_deform_enforce_weight, DivergenceFreeReg, DebugFolder, verbose);
    toc

    figure;imagescn(warped, [], [], [], 3)

    res_data = Matlab_gt_apply_deformation_field_reg_2D_series(single(data), dx, dy);
    res_mask = Matlab_gt_apply_deformation_field_reg_2D_series(single(mask), dx, dy);
    res_cine = Matlab_gt_apply_deformation_field_reg_2D_series(single(cine), dx, dy);
    figure; imagescn(data1-data2, [], [], [12]);
    figure; imagescn(data1-res_data, [], [], [12]);
    figure; imagescn([dx, dy]);

    if (numel(size(mask)) == 2)
        [rad_strain, circ_strain] = strain_Fn(mask, dx, dy);
    else
        [rad_strain, circ_strain] = strain_Fn(mask(:, :, phase), dx, dy);
    end
    
    % rad_strain(:, :, 1) = cine(:, :, 1);
    % radial_strain(:, :, i, :) = rad_strain;
    % circumfrential_strain(:, :, i, :) = circ_strain; 

end

end