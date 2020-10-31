function prepare_ai_cmr_training_LGE(resDir, aiDir, pt_ids, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_num, check_processed, plot_flag, local_pic_dir)
% prepare_ai_cmr_training_LGE(resDir, aiDir, files_record_picked)

mkdir(aiDir)

if(nargin>=11)
    pic_dir = fullfile(aiDir, local_pic_dir);
else    
    pic_dir = fullfile(aiDir, 'jpg_pics');
end

mkdir(pic_dir)
mkdir(fullfile(pic_dir, 'ch4'))
mkdir(fullfile(pic_dir, 'ch2'))
mkdir(fullfile(pic_dir, 'ch3'))
mkdir(fullfile(pic_dir, 'sax'))
mkdir(fullfile(pic_dir, 'ch4_numpy'))
mkdir(fullfile(pic_dir, 'ch2_numpy'))
mkdir(fullfile(pic_dir, 'ch3_numpy'))
mkdir(fullfile(pic_dir, 'sax_numpy'))

visible_status = 'on';
if(~plot_flag)
    visible_status = 'off';
end

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    
    case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '4ch');
    case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '2ch');
    case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '3ch');
    case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'sa');

    if(numel(case_4ch)==0 || numel(case_2ch)==0 || numel(case_sax)==0)
        continue;
    end
        
    case_prefix = [case_4ch.study_dates(end,:) '_' pt_id];
    dst_dir = fullfile(aiDir, case_4ch.study_dates(end,:), pt_id);
    
    ch4_processed = 0;
    ch2_processed = 0;
    ch3_processed = 0;
    sax_processed = 0;
    
    if(check_processed)
        dst_dir_sax = fullfile(dst_dir, 'sax');
        if(numel(case_3ch)==0)
            if(exist(fullfile(dst_dir_sax, 'record_header.mat')) ...
                & exist(fullfile(pic_dir, 'ch2', [case_prefix '_ch2_' '1.jpg'])) ...
                & exist(fullfile(pic_dir, 'ch4', [case_prefix '_ch4_' '1.jpg'])) ...
                & exist(fullfile(pic_dir, 'sax', [case_prefix '_sax_' '6.jpg'])) ...
                & exist(fullfile(pic_dir, 'sax_original_numpy', [case_prefix '_sax_' '5.npy'])) )
            
                disp(['already processed ' pt_id]);

                continue;
            end
        else
            if(exist(fullfile(dst_dir_sax, 'record_header.mat')) ...
                    & exist(fullfile(pic_dir, 'ch2', [case_prefix '_ch2_' '1.jpg'])) ...
                    & exist(fullfile(pic_dir, 'ch3', [case_prefix '_ch3_' '1.jpg'])) ...
                    & exist(fullfile(pic_dir, 'ch4', [case_prefix '_ch4_' '1.jpg'])) ...
                    & exist(fullfile(pic_dir, 'sax', [case_prefix '_sax_' '5.jpg'])) ...
                    & exist(fullfile(pic_dir, 'sax_original_numpy', [case_prefix '_sax_' '6.npy'])) )

                disp(['already processed ' pt_id]);

                continue;
            end
        end
        
        if(exist(fullfile(pic_dir, 'ch4_original_numpy', [case_prefix '_ch4_' '1.npy'])))
            ch4_processed = 1;
        end
        
        if(exist(fullfile(pic_dir, 'ch2_original_numpy', [case_prefix '_ch2_' '1.npy'])))
            ch2_processed = 1;
        end
        
        if(exist(fullfile(pic_dir, 'ch3_original_numpy', [case_prefix '_ch3_' '1.npy'])))
            ch3_processed = 1;
        end
        
        if(exist(fullfile(pic_dir, 'sax_original_numpy', [case_prefix '_sax_' '6.npy'])))
            sax_processed = 1;
        end
    end
    
    try
        SLC = case_sax.headers{end}.encoding.encodingLimits.slice.maximum +1;
        sax_ind = size(case_sax, 1);
        if(SLC<6)
            SLC = case_sax.headers{1}.encoding.encodingLimits.slice.maximum +1;            
            sax_ind = 1;
        end

        % -----------------------------
        
        if(~ch4_processed)
            case_4ch_file_name = case_4ch.file_names{end};
            case_4ch_dir = fullfile(resDir, case_4ch.study_dates(end,:), case_4ch.file_names{end})
            [gt_4ch, gt_h_4ch, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_4ch_dir, series_num, 1);
            gt_4ch = squeeze(gt_4ch);
        end
        
        % -----------------------------
        
        if(~ch3_processed)
            if(numel(case_3ch)>0)
                case_3ch_file_name = case_3ch.file_names{end};
                case_3ch_dir = fullfile(resDir, case_3ch.study_dates(end,:), case_3ch.file_names{end})
                [gt_3ch, gt_h_3ch, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_3ch_dir, series_num, 1);
                gt_3ch = squeeze(gt_3ch);
            end
        end
        % -----------------------------
        
        if(~ch2_processed)
            case_2ch_file_name = case_2ch.file_names{end};
            case_2ch_dir = fullfile(resDir, case_2ch.study_dates(end,:), case_2ch.file_names{end})       
            [gt_2ch, gt_h_2ch, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_2ch_dir, series_num, 1);
            gt_2ch = squeeze(gt_2ch);
        end
        
        if(~sax_processed)
            case_sax_file_name = case_sax.file_names{sax_ind};
            case_sax_dir = fullfile(resDir, case_sax.study_dates(sax_ind,:), case_sax.file_names{sax_ind})       
            multi_slice_mode = 1;
            try   
                [gt_sax_all, gt_h_sax_all, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_sax_dir, series_num, 1);
                gt_sax_all = squeeze(gt_sax_all);
                multi_slice_mode = 0;
            catch
                rmdir(dst_dir, 's');
                continue;
            end
        end
        
        if(plot_flag)
            h1 = figure; imagescn(permute(gt_4ch, [2 1]), [], [], [8]);
            h2 = figure; imagescn(permute(gt_2ch, [2 1]), [], [], [8]);
            h3 = figure; imagescn(permute(gt_3ch, [2 1]), [], [], [8]);            
            h4 = figure; imagescn(permute(gt_sax_all, [2 1 3]), [], [], [8]);
        end
    catch
        continue;
    end
       
    try
        mkdir(dst_dir)

        % 4ch
        if(~ch4_processed)           
            prepare_one_view(case_prefix, dst_dir, pic_dir, 'ch4', gt_4ch, gt_h_4ch, NN_RO, NN_E1, dst_pixel_spacing, case_4ch_file_name, plot_flag, visible_status);
        end
        
        % 2ch
        if(~ch2_processed)
            prepare_one_view(case_prefix, dst_dir, pic_dir, 'ch2', gt_2ch, gt_h_2ch, NN_RO, NN_E1, dst_pixel_spacing, case_2ch_file_name, plot_flag, visible_status);
        end
        
        % 3ch
        if(~ch3_processed && numel(case_3ch)>0)
            prepare_one_view(case_prefix, dst_dir, pic_dir, 'ch3', gt_3ch, gt_h_3ch, NN_RO, NN_E1, dst_pixel_spacing, case_3ch_file_name, plot_flag, visible_status);
        end
    catch
        rmdir(dst_dir, 's');
        continue;
    end
        
    %sax
    if(~sax_processed)        
        try
            prepare_one_view(case_prefix, dst_dir, pic_dir, 'sax', gt_sax_all, gt_h_sax_all, NN_RO, NN_E1, dst_pixel_spacing, case_sax_file_name, plot_flag, visible_status);
        catch
            rmdir(dst_dir, 's');
            continue;
        end   
    end
    
    closeall
    closeall    
end

end

function [gt_view_sorted, gt_h_view_sorted] = sort_basal_to_apical(gt_view, gt_h_view)
    SLC = size(gt_view, 3);
    if(SLC==1)
        gt_view_sorted = gt_view;
        gt_h_view_sorted = gt_h_view;
        return
    end
    
    SLC_LOC = zeros(SLC, 1);
    for slc=1:SLC
        SLC_LOC(slc) = gt_h_view(slc).slice_location;
    end
    
    [SL, ind] = sort(SLC_LOC);
    gt_view_sorted = gt_view(:,:,ind);
    gt_h_view_sorted = gt_h_view(ind,1);
end

function prepare_one_view(case_prefix, dst_dir, pic_dir, view_str, gt_view, gt_h_view, NN_RO, NN_E1, dst_pixel_spacing, file_name, plot_flag, visible_status)
    dst_dir_view = fullfile(dst_dir, view_str);
    mkdir(dst_dir_view);

    [gt_view, gt_h_view] = sort_basal_to_apical(gt_view, gt_h_view);
    
    RO = size(gt_view,1)
    E1 = size(gt_view,2)
    ps = max(gt_h_view(1,1).FOV)/max(RO,E1)

    WC = gt_h_view(1,1).window_center;
    WW = gt_h_view(1,1).window_width;
    
    new_RO = round(RO*ps/dst_pixel_spacing(1));
    new_E1 = round(E1*ps/dst_pixel_spacing(2));

    gt_view(find(isnan(gt_view))) = 0;
    data = Matlab_gt_resize_2D_image(double(gt_view), new_RO, new_E1, 5);
    if(plot_flag)
        figure; imagescn(permute(data, [2 1 3]), [WC-WW/2 WC+WW/2], [], [8]);
    end

    gt_view = permute(gt_view, [2 1 3]);
    data = permute(data, [2 1 3]);
    size(data)

    RO = size(data,1);
    E1 = size(data,2);
    PHS = size(data,3);

    data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), PHS);

    RO = size(data,1);
    E1 = size(data,2);
    SLC = size(data,3);

    data_resized_training = data;
    h2 = figure('visible', visible_status); imagescn(data_resized_training, [WC-WW/2 WC+WW/2], [], [8]);
    set(h2, 'visible', 'off')
    header = CreateGtImageHeader(data_resized_training);
    Matlab_gt_write_analyze(single(data_resized_training), header, fullfile(dst_dir_view, 'data_resized_training'));
   
    writeNPY(single(gt_view), fullfile(dst_dir_view, 'data_acq.npy'));
    writeNPY(single(data_resized_training), fullfile(dst_dir_view, 'data_resized_training.npy'));
    writeNPY(single(data), fullfile(dst_dir_view, 'data.npy'));
    saveas(h2, fullfile(dst_dir_view, 'data_resized_training'), 'jpg');

    pixel_spacing = ps;
    save(fullfile(dst_dir_view, 'record'), 'gt_h_view', 'dst_pixel_spacing', 'pixel_spacing', 'file_name');
    save(fullfile(dst_dir_view, 'record_header'), 'gt_h_view', 'dst_pixel_spacing', 'pixel_spacing', 'file_name');

    for slc=1:SLC
        h2 = figure('visible', visible_status); 
        
        WC = gt_h_view(slc).window_center;
        WW = gt_h_view(slc).window_width;
    
        imagescn(data_resized_training(:,:,slc), [WC-WW/2 WC+WW/2], [], [8]);
        set(h2, 'visible', 'off')
        saveas(h2, fullfile(dst_dir_view, [case_prefix '_' view_str '_' num2str(slc) '.jpg']), 'jpg');
        saveas(h2, fullfile(pic_dir, view_str, [case_prefix '_' view_str '_' num2str(slc) '.jpg']), 'jpg');    
        writeNPY(single(data_resized_training(:,:,slc)), fullfile(pic_dir, [view_str '_numpy'], [case_prefix '_' view_str '_' num2str(slc) '.npy']));
    end
end