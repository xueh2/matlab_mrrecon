function prepare_ai_cmr_training_MOLLI(resDir, aiDir, pt_ids, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_num, series_num_map, check_processed, plot_flag, local_pic_dir)
% prepare_ai_cmr_training_MOLLI(resDir, aiDir, files_record_picked)

mkdir(aiDir)

if(nargin>=11)
    pic_dir = fullfile(aiDir, local_pic_dir);
else    
    pic_dir = fullfile(aiDir, 'jpg_pics');
end

mkdir(pic_dir)
mkdir(fullfile(pic_dir, 'pre'))
mkdir(fullfile(pic_dir, 'post'))
mkdir(fullfile(pic_dir, 'pre_numpy'))
mkdir(fullfile(pic_dir, 'post_numpy'))
mkdir(fullfile(pic_dir, 'pre_map'))
mkdir(fullfile(pic_dir, 'post_map'))
mkdir(fullfile(pic_dir, 'pre_map_numpy'))
mkdir(fullfile(pic_dir, 'post_map_numpy'))
mkdir(fullfile(pic_dir, 'post_image'))
mkdir(fullfile(pic_dir, 'post_image_numpy'))

visible_status = 'on';
if(~plot_flag)
    visible_status = 'off';
end

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    
    case_pre = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'pre');
    case_post = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'post');

    if(numel(case_pre)==0 || numel(case_post)==0)
        continue;
    end
        
    case_prefix = [case_pre.study_dates(end,:) '_' pt_id];
    dst_dir = fullfile(aiDir, case_pre.study_dates(end,:), pt_id);
       
    for v=1:2        
        if(v==1)
            case_used = case_pre;
            pic_dir_used = fullfile(pic_dir, 'pre');
        else
            case_used = case_post;
            pic_dir_used = fullfile(pic_dir, 'post');
        end
        
        for k=1:size(case_used, 1)
            dst_dir_case = fullfile(dst_dir, case_used.file_names{k});
            if(check_processed ... 
                & exist(fullfile(dst_dir_case, 'record_header.mat')) ...
                & exist(fullfile(pic_dir_used, [case_used.file_names{k} '_1.jpg'])) ...
                & exist(fullfile([pic_dir_used '_numpy'], [case_used.file_names{k} '_1.npy'])) )

                disp(['already processed ' pt_id]);

                continue;
            end

            % process this case
            try
                case_dir = fullfile(resDir, case_used.study_dates(k,:), case_used.file_names{k})
                [gt_moco, gt_h_moco, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, series_num, 1);
                gt_moco = squeeze(gt_moco);

                [gt_map, gt_h_map, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, series_num_map, 1);
                gt_map = squeeze(gt_map);
            catch
                try
                    case_dir = fullfile(resDir, case_used.study_dates(k,:), case_used.file_names{k})
                    [gt_moco, gt_h_moco, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, series_num-1, 1);
                    gt_moco = squeeze(gt_moco);

                    [gt_map, gt_h_map, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, series_num_map-1, 1);
                    gt_map = squeeze(gt_map);
                catch
                    continue;
                end
            end
            
            I = [];
            for slc=1:size(gt_moco, 3)            
                im = gt_moco(:,:,slc);
                im = flipdim(flipdim(im, 2), 1);
                if(slc==1)
                    ImageOrientationPatient = [gt_h_moco(1).read_dir gt_h_moco(1).phase_dir gt_h_moco(1).slice_dir];
                    [rotate90, flip] = NormOrientation_rev3 (ImageOrientationPatient);
                end
                [im_cmr_view, gt_h] = apply_normorientation(im, rotate90, flip, gt_h_moco(slc));
    %             figure; imagescn(im_cmr_view);
                I(:,:,slc) = im_cmr_view;
                gt_h_moco(slc) = gt_h;
            end
            gt_moco = I;

            im = gt_map;
            im = flipdim(flipdim(im, 2), 1);
            ImageOrientationPatient = [gt_h_map(1).read_dir gt_h_map(1).phase_dir gt_h_map(1).slice_dir];
            [rotate90, flip] = NormOrientation_rev3 (ImageOrientationPatient);
            [im_cmr_view, gt_h] = apply_normorientation(im, rotate90, flip, gt_h_map(1));
    %         figure; imagescn(im_cmr_view);
            gt_map = im_cmr_view;
            gt_h_map(1) = gt_h;

%             prepare_one_view(case_used.file_names{k}, dst_dir, pic_dir_used, [], gt_moco, gt_h_moco, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, v==1);
            prepare_one_view(case_used.file_names{k}, dst_dir, [pic_dir_used '_map'], 't1_map', gt_map, gt_h_map, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, v==2);
            
            if(v==2)
                set_used = 5;
                if(set_used>size(gt_moco, 3))
                    set_used = size(gt_moco, 3);
                end
                prepare_one_view(case_used.file_names{k}, dst_dir, [pic_dir_used '_image'], 't1_post_image', gt_moco(:,:,set_used), gt_h_moco(set_used), NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, v==2);
            end
        end
    end
    
    closeall
    closeall    
end

end

function prepare_one_view(case_prefix, dst_dir, pic_dir, view_str, gt_view, gt_h_view, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, is_post)
    if(isempty(view_str))
        dst_dir_view = fullfile(dst_dir, case_prefix);
    else
        dst_dir_view = fullfile(dst_dir, [case_prefix '__' view_str]);
    end
    mkdir(dst_dir_view);
   
    RO = size(gt_view,1);
    E1 = size(gt_view,2);
    ps = max(gt_h_view(1,1).FOV)/max(RO,E1);

    WC = gt_h_view(1,1).window_center;
    WW = gt_h_view(1,1).window_width;

    if(WC<0)
        WC = median(gt_view(:));
        WW =  0.25 * mean(gt_view(:));
    end
    
    if(~isempty(view_str) & strcmp(view_str,'t1_map') & is_post)
        WC = 550;
        WW = 650;
    end
    
    new_RO = round(RO*ps/dst_pixel_spacing(1));
    new_E1 = round(E1*ps/dst_pixel_spacing(2));

    gt_view(find(isnan(gt_view))) = 0;
    data = Matlab_gt_resize_2D_image(double(gt_view), new_RO, new_E1, 5);
    if(plot_flag)
        figure; imagescn(data, [WC-WW/2 WC+WW/2], [], [8]);
    end

    RO = size(data,1);
    E1 = size(data,2);
    SET = size(data,3);

    data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), SET);

    RO = size(data,1);
    E1 = size(data,2);
    SLC = size(data,3);

    data_resized_training = data;
    h2 = figure('visible', visible_status); imagescn(data_resized_training, [WC-WW/2 WC+WW/2], [], [8]);
    if(~isempty(view_str) & strcmp(view_str,'t1_map'))
        p = T1ColorMap;
        colormap(p);
    end
    set(h2, 'visible', 'off')
    header = CreateGtImageHeader(data_resized_training);
    Matlab_gt_write_analyze(single(data_resized_training), header, fullfile(dst_dir_view, 'data_resized_training'));

    writeNPY(single(gt_view), fullfile(dst_dir_view, 'data_acq.npy'));
    writeNPY(single(data_resized_training), fullfile(dst_dir_view, 'data_resized_training.npy'));
    writeNPY(single(data), fullfile(dst_dir_view, 'data.npy'));
    saveas(h2, fullfile(dst_dir_view, 'data_resized_training'), 'jpg');

    pixel_spacing = ps;
    save(fullfile(dst_dir_view, 'record'), 'gt_h_view', 'dst_pixel_spacing', 'pixel_spacing');
    save(fullfile(dst_dir_view, 'record_header'), 'gt_h_view', 'dst_pixel_spacing', 'pixel_spacing');
    
    for s=1:SET
        
%         if(~isempty(view_str) & strcmp(view_str,'t1_post_image') & is_post)
%             pp = data_resized_training(:,:,s);
%             max_v = max(pp(:));
%             a2 = adapthisteq(uint8(255 * pp/max_v));
%             h2 = figure('visible', visible_status); imagescn(a2, [0*median(a2(:)) 1.5*median(a2(:))]);
%         else
            h2 = figure('visible', visible_status); imagescn(data_resized_training(:,:,s), [WC-WW/2 WC+WW/2], [], [8]);
            if(~isempty(view_str) & strcmp(view_str,'t1_map'))
                p = T1ColorMap;
                colormap(p);
            end
%         end
        saveas(h2, fullfile(pic_dir, [case_prefix '_' num2str(s) '.jpg']), 'jpg');
        writeNPY(single(data_resized_training(:,:,s)), fullfile([pic_dir '_numpy'], [case_prefix '_' num2str(s) '.npy']));
    end
    closeall
end