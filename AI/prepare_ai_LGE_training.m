function prepare_ai_LGE_training(resDir, aiDir, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_num, permute_data)
% prepare_ai_LGE_training(resDir, aiDir, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_num, permute_data)

mkdir(aiDir)

pt_ids = unique(files_record_picked(:, 3));

for pt=1:size(pt_ids, 1)
    
    closeall
    
    pt_id = pt_ids.patientIDs{pt};
    
    disp([num2str(pt) ' out of ' num2str(size(pt_ids,1)) ' - ' pt_id]);
    
    case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '4ch');
    case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '2ch');
    case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'sa');

    if(numel(case_4ch)==0 || numel(case_2ch)==0 || numel(case_sax)==0)
        continue;
    end
        
    dst_dir = fullfile(aiDir, case_4ch.study_dates(end,:), pt_id);
    if(exist(dst_dir) & exist(fullfile(dst_dir, 'sax', 'Cine_resized_training.jpg')))
        continue;
    end
    
    try
        SLC = case_sax.headers{end}.encoding.encodingLimits.slice.maximum +1;
        sax_ind = size(case_sax, 1);
        if(SLC<6)
            SLC = case_sax.headers{1}.encoding.encodingLimits.slice.maximum +1;            
            sax_ind = 1;
        end

        case_4ch_dir = fullfile(resDir, case_4ch.study_dates(end,:), case_4ch.file_names{end})
        [gt_4ch, gt_h_4ch, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_4ch_dir, series_num, 1);
        gt_4ch = squeeze(gt_4ch);
        
        % -----------------------------
        
        case_2ch_dir = fullfile(resDir, case_2ch.study_dates(end,:), case_2ch.file_names{end})       
        [gt_2ch, gt_h_2ch, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_2ch_dir, series_num, 1);
        gt_2ch = squeeze(gt_2ch);

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

        if(size(gt_sax_all,4)~= size(gt_2ch,3))
            multi_slice_mode = 1;
        end

        if(numel(size(gt_sax_all))==4)
            if(size(gt_sax_all,4)~= size(gt_2ch,3))
                rmdir(dst_dir, 's');
                continue;
            end
        end
        
        if(permute_data)
            h1 = figure; imagescn(permute(gt_4ch, [2 1 3]), [], [], [8]);
            h2 = figure; imagescn(permute(gt_2ch, [2 1 3]), [], [], [8]);
            h3 = figure; imagescn(permute(gt_sax_all, [2 1 3]), [], [], [8]);
        else
            h1 = figure; imagescn(gt_4ch, [], [], [8]);
            h2 = figure; imagescn(gt_2ch, [], [], [8]);
            h3 = figure; imagescn(gt_sax_all, [], [], [8]);
        end
    catch
        continue;
    end
    
    mkdir(dst_dir)

    % 4ch
    dst_dir_4ch = fullfile(dst_dir, 'ch4');
    mkdir(dst_dir_4ch);
           
    RO = size(gt_4ch,1)
    E1 = size(gt_4ch,2)
    ps = max(gt_h_4ch(1,1,1).FOV)/max(RO,E1)
    
    new_RO = round(RO*ps/dst_pixel_spacing(1));
    new_E1 = round(E1*ps/dst_pixel_spacing(2));
    
    data = Matlab_gt_resize_2D_image(double(gt_4ch), new_RO, new_E1, 5);   
    
    if(permute_data)
        data = permute(data, [2 1 3]);
    end
    size(data)
   
    RO = size(data,1);
    E1 = size(data,2);
    PHS = size(data,3);
    
    data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), PHS);
   
    RO = size(data,1);
    E1 = size(data,2);
    PHS = size(data,3);
    
    [ro_s, ro_e, e1_s, e1_e] = detect_heart_region_LGE_roi(data, NN_RO, NN_E1);
    
    Cine_resized_training = data(ro_s:ro_e, e1_s:e1_e, :, :);
    h2 = figure; imagescn(Cine_resized_training, [], [], [8]);
    
    header = CreateGtImageHeader(Cine_resized_training);
    Matlab_gt_write_analyze(single(Cine_resized_training), header, fullfile(dst_dir_4ch, 'Cine_resized_training'));
   
    writeNPY(single(gt_4ch), fullfile(dst_dir_4ch, 'data_acq.npy'));   
    writeNPY(single(Cine_resized_training), fullfile(dst_dir_4ch, 'Cine_resized_training.npy'));
    
    writeNPY(single(data), fullfile(dst_dir_4ch, 'data.npy'));
    writeNPY(single([ro_s, ro_e, e1_s, e1_e]), fullfile(dst_dir_4ch, 'roi.npy'));
    
    saveas(h2, fullfile(dst_dir_4ch, 'Cine_resized_training'), 'jpg');
    
    pixel_spacing = ps;
    save(fullfile(dst_dir_4ch, 'record'), 'gt_4ch', 'gt_h_4ch', 'dst_pixel_spacing', 'pixel_spacing');
    
    % 2ch
    dst_dir_2ch = fullfile(dst_dir, 'ch2');
    mkdir(dst_dir_2ch);
       
    RO = size(gt_2ch,1)
    E1 = size(gt_2ch,2)
    ps = max(gt_h_2ch(1,1,1).FOV)/max(RO,E1)
    
    new_RO = round(RO*ps/dst_pixel_spacing(1));
    new_E1 = round(E1*ps/dst_pixel_spacing(2));
    
    data = Matlab_gt_resize_2D_image(double(gt_2ch), new_RO, new_E1, 5);
    
    if(permute_data)
        data = permute(data, [2 1 3]);
    end
    size(data)
   
    RO = size(data,1);
    E1 = size(data,2);
    PHS = size(data,3);
    
    data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), PHS);
   
    RO = size(data,1);
    E1 = size(data,2);
    PHS = size(data,3);
    
    [ro_s, ro_e, e1_s, e1_e] = detect_heart_region_LGE_roi(data, NN_RO, NN_E1);
    
    Cine_resized_training = data(ro_s:ro_e, e1_s:e1_e, :, :);
    h2 = figure; imagescn(Cine_resized_training, [], [], [8], 3);
    
    header = CreateGtImageHeader(Cine_resized_training);
    Matlab_gt_write_analyze(single(Cine_resized_training), header, fullfile(dst_dir_2ch, 'Cine_resized_training'));
   
    writeNPY(single(Cine_resized_training), fullfile(dst_dir_2ch, 'Cine_resized_training.npy'));   
    writeNPY(single(gt_2ch), fullfile(dst_dir_2ch, 'data_acq.npy'));   
    writeNPY(single(data), fullfile(dst_dir_2ch, 'data.npy'));
    writeNPY(single([ro_s, ro_e, e1_s, e1_e]), fullfile(dst_dir_2ch, 'roi.npy'));
    
    saveas(h2, fullfile(dst_dir_2ch, 'Cine_resized_training'), 'jpg');
    
    pixel_spacing = ps;
    save(fullfile(dst_dir_2ch, 'record'), 'gt_2ch', 'gt_h_2ch', 'dst_pixel_spacing', 'pixel_spacing');
    
    %sax
    dst_dir_sax = fullfile(dst_dir, 'sax');
    mkdir(dst_dir_sax);
            
    % resample sax
    RO = size(gt_sax_all,1)
    E1 = size(gt_sax_all,2)
    ps = max(gt_h_sax_all(1,1,1).FOV)/max(RO,E1)
    
    new_RO = round(RO*ps/dst_pixel_spacing(1));
    new_E1 = round(E1*ps/dst_pixel_spacing(2));
    
    data = Matlab_gt_resize_2D_image(double(gt_sax_all), new_RO, new_E1, 5);
    
    if(permute_data)
        data = permute(data, [2 1 4 3]);
    end
    size(data)
   
    RO = size(data,1);
    E1 = size(data,2);
    PHS = size(data,3);
    
    data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), PHS);
   
    RO = size(data,1);
    E1 = size(data,2);
    PHS = size(data,3);
    
    [ro_s, ro_e, e1_s, e1_e] = detect_heart_region_LGE_roi(data, NN_RO, NN_E1);
        
    Cine_resized_training = data(ro_s:ro_e, e1_s:e1_e, :, :);
    h2 = figure; imagescn(Cine_resized_training, [], [], [8]);
    
    header = CreateGtImageHeader(Cine_resized_training);
    Matlab_gt_write_analyze(single(Cine_resized_training), header, fullfile(dst_dir_sax, 'Cine_resized_training'));
   
    writeNPY(single(Cine_resized_training), fullfile(dst_dir_sax, 'Cine_resized_training.npy'));   
    writeNPY(single(gt_sax_all), fullfile(dst_dir_sax, 'data_acq.npy'));    
    writeNPY(single(data), fullfile(dst_dir_sax, 'data.npy'));
    writeNPY(single([ro_s, ro_e, e1_s, e1_e]), fullfile(dst_dir_sax, 'roi.npy'));
    
    saveas(h2, fullfile(dst_dir_sax, 'Cine_resized_training'), 'jpg');
    
    pixel_spacing = ps;
    save(fullfile(dst_dir_sax, 'record'), 'gt_sax_all', 'gt_h_sax_all', 'dst_pixel_spacing', 'pixel_spacing');
    
    closeall
    closeall    
end
