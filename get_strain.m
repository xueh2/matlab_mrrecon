%% load mask from file
clear all
close all
cd C:\gtuser\strain_2
filename = '229201050_229201055_362_20190718'
load([filename, '.mat'])
load([filename, '_mask.mat'])

%%

cine = double(rt_cine);
mask = rt_mask;
PHS = size(cine,3);

[rad_strain, circ_strain] = strainWrapper(cine, mask, @numericalStrain);

figure
if (numel(size(mask)) == 2)
    imagescn([rad_strain(:, :, 2:PHS).*repmat(mask, 1, 1, PHS-1), circ_strain(:, :, 2:PHS).*repmat(mask, 1, 1, PHS-1)], [-0.1, 0.5], [], [], 3)
else
    imagescn([rad_strain(:, :, 2:PHS).*mask(:, :, 2:PHS), circ_strain(:, :, 2:PHS).*mask(:, :, 2:PHS)], [-0.1, 0.5], [], [], 3)
end
% save([filename, 'rt_strain.mat'], 'rad_strain', 'circ_strain');


%% get centered masks
clear all
close all
cd C:\gtuser\
filename = 'Tag4-6-2'
cd(filename)

load 229201082_229201087_478_20190718

mask = rt_cine_epi_mask - rt_cine_endo_mask;
cine = rt_cine;
figure; imagescn(cine, [], [], [], 3);

%% 

cine = double(cine);
mask = double(mask);
mask(mask < 0) = 0;
mask(isnan(mask)) = 0;
cine(isnan(cine)) = 0;
cine(isnan(mask)) = 0;

figure; imagescn(cine, [], [], [], 3)
figure; imagescn(mask, [], [], [], 3)

% deformation_data = zeros(mask_size);
% for i = 1:mask_size(3)
%     mask_stats = regionprops(mask(:, :, i), 'centroid');
%     [X, Y] = meshgrid(1:mask_size(1), 1:mask_size(2));
%     dists = sqrt((X - mask_stats.Centroid(1)).^2 + (Y - mask_stats.Centroid(2)).^2);
%     mask_dists = 5*mask(:, :, i).*dists;
%     mask_dists = mask_dists/max(max(max(mask_dists)));
%     mask_dists(mask_dists > 0) = mask_dists(mask_dists > 0) - min(min((mask_dists(mask_dists > 0))/max(max(mask_dists(mask_dists> 0)))));
%     inv_mask_dists = mask(:, :, i).*(1 - mask_dists);
%     inv_mask_dists = inv_mask_dists/max(max(max(inv_mask_dists)));
%     inv_mask_dists(inv_mask_dists > 0) = inv_mask_dists(inv_mask_dists > 0) - min(min((inv_mask_dists(inv_mask_dists > 0))/max(max(inv_mask_dists(inv_mask_dists> 0)))));
%     deformation_data(:, :, i) = (((13*mask_dists).^0.5 + (13*inv_mask_dists).^0.5) + 1);
% end
% deformation_data = deformation_data.*cine;
% 
% % gaussian mask
% mask_neg = mask;
% mask_neg(mask == 0) = -1;
% deformation_data = imgaussfilt(mask_neg,2);
% deformation_data = deformation_data - min(min(min(deformation_data(mask == 1))));
% deformation_data(mask == 0) = 0;
% imagescn(deformation_data, [], [], [], 3)
% deformation_data = (1+deformation_data * 5) .* cine;
% imagescn(deformation_data, [], [], [], 3)
% 
% 
% % inv gaussian mask
% mask_neg = mask;
% mask_neg(mask == 0) = -1;
% deformation_data = 1/(imgaussfilt(mask_neg,2)+0.2);
% deformation_data(isnan(deformation_data)) = 1;
% deformation_data(deformation_data > 5) = 1.5;
% deformation_data(deformation_data < 0) = 0;
% deformation_data(mask == 0) = 0;
% figure;imagescn(deformation_data, [0, 4], [], [], 3)
% deformation_data = (1+deformation_data*4) .* cine;
% figure;imagescn(deformation_data, [], [], [], 3)
% 
% deformation_data = mask*200+cine;
% deformation_data = mask;

deformation_data = rt_cine + 1000 * mask;
shift_centroid = true;
add_centroid = true;
imagescn(deformation_data, [], [], [], 3);
if shift_centroid
    [deformation_data, shifted_cine, shifted_mask, old_deformation_data] = fixCentroid(deformation_data, cine, mask, add_centroid);
end
figure; imagescn([old_deformation_data], [], [], [], 3)
figure; imagescn([deformation_data], [], [], [], 3)
[rad_strain, circ_strain, res_data] =  strainWrapper(deformation_data, shifted_cine, shifted_mask, @numericalStrain);
% [rad_strain, circ_strain] = numericalStrain(mask(:, :, 10), dx, dy);
figure;imagescn([rad_strain.*shifted_mask(:, :, :), circ_strain.*shifted_mask(:, :, :)], [-0.8, 0.8], [], [], 3)
figure;imagescn([rad_strain, circ_strain], [-0.8, 0.8], [], [], 3)

%% load all data and get results from all
clear all
closeall
files = {'C:\gtuser\20190716_Tag1__study_\', 'C:\gtuser\20190717_Tag2__study_\', 'C:\gtuser\20190717_Tag3__study_\'};
rt_tagged_bh = {'00024', '00027', '00112'; '00016', '00019', '00112'};


%% get rt, tagged, bh data

for set = 1:3
    f = files{1, set};
    cd(f);
    
    for i = 1:3
        folder = rt_tagged_bh{set, i};
        cd(folder);
        phs = numel(dir('*.dcm'));
        
        for p = 1:phs
            if (set == 2) && (i == 3)
                temp = dicomread([num2str(p, '%04u'), '_2.dcm']);
            else
                temp = dicomread([num2str(p, '%04u'), '.dcm']);
            end
            
            if p == 1
                im_sz = size(temp);
                im = zeros(im_sz(1), im_sz(2), phs);
            end
            im(:, :, p) = temp;
        end

        if i == 1
            rt = im;
        elseif i == 2
            tagged = im;
        elseif i == 3
            bh = im;
        end
        cd ../
    end
    
    save('data_all.mat', 'rt', 'tagged', 'bh')
    cd ../
end

%% get contours for all three
file_idx = 2;
f = files{1, file_idx}

cd(f)
load('data_all.mat')
figure; imagescn(rt, [], [], [], 3); title('bh'); pause
figure; imagescn(tagged, [], [], [], 3); title('tagged'); pause
figure; imagescn(bh, [], [], [], 3); title('rt'); pause

%% get masks
for set = 2:2
    f = files{1, set};
    cd(f);
    load('data_all.mat')
    labels = {'rt', 'tagged', 'bh'};
    for i = 1:3        
        load(['roi_', labels{i}, '.mat'])
        epi_x = ROI_info_table(1).ROI_x_coordinates;
        endo_x = ROI_info_table(2).ROI_x_coordinates;
        epi_y = ROI_info_table(1).ROI_y_coordinates;
        endo_y = ROI_info_table(2).ROI_y_coordinates;
        
        if i == 1
            cine_sz = size(rt);
        elseif i == 2
            cine_sz = size(tagged);
        elseif i == 3
            cine_sz = size(bh);
        end
        
        endo = poly2mask(endo_x,endo_y,cine_sz(1), cine_sz(2));
        epi = poly2mask(epi_x,epi_y,cine_sz(1), cine_sz(2));

        mask = epi-endo;
        figure
        imagescn(mask)
        
        if i == 1
            rt_mask = mask;
        elseif i == 2
            tagged_mask = mask;
        elseif i == 3
            bh_mask = mask;
        end
    end
    
    save('masks_all.mat', 'rt_mask', 'tagged_mask', 'bh_mask')
    cd ../
end
sum(sum(sum(sum(isnan(masks_all)))))


%% get strain maps

for set = 1:1
    f = files{1, set};
    cd(f);
    load('data_all.mat')
    load('masks_all.mat')
    for i = 1:3
        
        if i == 1
            cine = double(rt);
            mask = rt_mask; 
        elseif i == 2
            cine = double(tag);
            mask = tagged_mask;
        elseif i == 3
            cine = double(bh);
            mask = bh_mask;
        end

        dissimilarity = 'LocalCCR';
        level = 4;
        max_iter_num_pyramid_level = [64 100 100 100];
        regularization_hilbert_strength = [16 16 16 16]/2;

        LocalCCR_sigmaArg = [2 2 2 2];
        BidirectionalReg = 1;
        dissimilarity_thres = 0;
        div_num = 5;
        inverse_deform_enforce_iter = 10;
        inverse_deform_enforce_weight = 0.5;
        DivergenceFreeReg = 0;
        DebugFolder = [];
        verbose = 0;

        key_frame = 0;

        data = cine;

        RO = size(data,1);
        E1 = size(data,2);
        PHS = size(data,3);

        data1 = data(:,:,1:end-1);
        data2 = data(:,:,2:end);

        % Get deformation field of im i to im i + 1

        tic
        [dx, dy, warped, dxInv, dyInv] = Matlab_gt_deformation_field_reg_2D_series_pairwise(single(data1), single(data2), dissimilarity, level, max_iter_num_pyramid_level, regularization_hilbert_strength, LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, div_num, inverse_deform_enforce_iter, inverse_deform_enforce_weight, DivergenceFreeReg, DebugFolder, verbose);
        toc

        % figure
        % imagescn([dx, dy], [], [], [], 3)

        res = Matlab_gt_apply_deformation_field_reg_2D_series(single(data2), dx, dy);

        dx1 = zeros(RO, E1, PHS);
        dx1(:,:,1) = 0;
        dx1(:,:,2:end) = dx;

        dy1 = zeros(RO, E1, PHS);
        dy1(:,:,1) = 0;
        dy1(:,:,2:end) = dy;

        figure; imagescn(data, [], [], [], 3)

        % figure
        % imagescn([dx1, dy1], [], [], [], 3)

        tic
        [dx_out, dy_out] = Matlab_gt_concatenate_deform_fields(single(dx1), single(dy1), key_frame); 
        toc

        res = zeros(size(data));
        [Y, X] = meshgrid(1:E1, 1:RO);

        for p = 1:PHS
            res(:, :, p) = interp2(data(:,:,p), dy_out(:,:,p) + Y, dx_out(:,:,p) + X);
        end
        figure; imagescn(res, [], [], [], 3)

        % figure
        % imagescn([dx_out, dy_out], [], [], [], 3)

        [rad_strain, circ_strain] = numericalStrain(mask, dx_out, dy_out);
        % rad_strain(:, :, 1) = cine(:, :, 1);
        % radial_strain(:, :, i, :) = rad_strain;
        % circumfrential_strain(:, :, i, :) = circ_strain; 


        figure
        imagescn(rad_strain(:, :, 2:PHS).*repmat(mask, 1, 1, PHS-1), [-0.1, 0.5], [], [], 3)
        figure
        imagescn(circ_strain(:, :, 2:PHS).*repmat(mask, 1, 1, PHS-1), [-0.1, 0.5], [], [], 3)

        if i == 1
            rad_strain_rt = rad_strain;
            circ_strain_rt = circ_strain;
            moco_rt = res;
        elseif i == 2
            rad_strain_tagged = rad_strain;
            circ_strain_tagged = circ_strain;
            moco_tagged = res;
        elseif i == 3
            rad_strain_bh = rad_strain;
            circ_strain_bh = circ_strain;
            moco_bh = res;
        end  
    end  
    save('radial_strain.mat', 'rad_strain_rt', 'rad_strain_tagged', 'rad_strain_bh');
    save('circ_strain.mat', 'circ_strain_rt', 'circ_strain_tagged', 'circ_strain_bh');
    save('moco.mat', 'moco_rt', 'moco_tagged', 'moco_bh')
end


