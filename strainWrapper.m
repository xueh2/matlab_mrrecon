function [rad_strain, circ_strain, res_data] = strainWrapper(data, cine, mask, strain_Fn)


dissimilarity = 'LocalCCR';
level = 4;
max_iter_num_pyramid_level = [200 200 200 200];
regularization_hilbert_strength = [16 16 16 16]/2;
LocalCCR_sigmaArg = [2 2 2 2];
BidirectionalReg = 0;
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

data1 = repmat(data(:,:,10), [1 1 PHS]);
data2 = data;
tic
[dx, dy, warped, dxInv, dyInv] = Matlab_gt_deformation_field_reg_2D_series_pairwise(double(data1), double(data2), dissimilarity, level, max_iter_num_pyramid_level, regularization_hilbert_strength, LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, div_num, inverse_deform_enforce_iter, inverse_deform_enforce_weight, DivergenceFreeReg, DebugFolder, verbose);
toc

% dx1 = zeros(RO, E1, PHS);
% dx1(:,:,1) = 0;
% dx1(:,:,2:end) = dx;
% 
% dy1 = zeros(RO, E1, PHS);
% dy1(:,:,1) = 0;
% dy1(:,:,2:end) = dy;

% figure
% imagescn([dx1, dy1], [], [], [], 3)

% tic
% [dx_out, dy_out] = Matlab_gt_concatenate_deform_fields(double(dx1), double(dy1), key_frame); 
% toc

res_data = Matlab_gt_apply_deformation_field_reg_2D_series(single(data), dx, dy);
res_mask = Matlab_gt_apply_deformation_field_reg_2D_series(single(mask), dx, dy);
res_cine = Matlab_gt_apply_deformation_field_reg_2D_series(single(cine), dx, dy);
figure; imagescn(data1-data2, [], [], [12]);
figure; imagescn(data1-res_data, [], [], [12]);
figure; imagescn([dx, dy]);
% 
% res_data = Matlab_gt_apply_deformation_field_reg_2D_series(single(data), dx, dy);
% res_mask = Matlab_gt_apply_deformation_field_reg_2D_series(single(mask), dx, dy);
% res_cine = Matlab_gt_apply_deformation_field_reg_2D_series(single(cine), dx, dy);

[rad_strain, circ_strain] = strain_Fn(mask(:, :, 10), dx, dy);
% rad_strain(:, :, 1) = cine(:, :, 1);
% radial_strain(:, :, i, :) = rad_strain;
% circumfrential_strain(:, :, i, :) = circ_strain; 

end