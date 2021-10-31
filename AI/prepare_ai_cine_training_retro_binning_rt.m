function prepare_ai_cine_training_retro_binning_rt(dataDir, resDir, aiDir, pt_ids, files_record_picked, dst_pixel_spacing, series_num, gmap_num, record, check_processed, plot_flag, local_pic_dir)
% prepare_ai_cine_training_retro_binning_rt(dataDir, resDir, aiDir, pt_ids, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_num, gmap_num, record, check_processed, plot_flag, local_pic_dir)

mkdir(aiDir)

if(~isempty(local_pic_dir))
    pic_dir = fullfile(aiDir, local_pic_dir);
else    
    pic_dir = fullfile(aiDir, 'jpg_pics');
end

mkdir(pic_dir)

view_strs = {'ch4', 'ch2', 'ch3', 'sax', 'lvot', 'aov', 'rv', 'aorticarch'};

for k=1:numel(view_strs)
    view_str = view_strs{k};
    mkdir(fullfile(pic_dir, view_str));
    mkdir(fullfile(pic_dir, [view_str '_original']));
    mkdir(fullfile(pic_dir, [view_str '_numpy']));
    mkdir(fullfile(pic_dir, [view_str '_original_numpy']));
end

visible_status = 'on';
if(~plot_flag)
    visible_status = 'off';
end

% pt_ids = unique(files_record_picked(:, 3));

ED = 1;
ES = 13;

useMask = 1;
boxFilterSize = 25;
noisebackground = 30;
thresRatioForNoise = 10;
    
for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};       
    
    if(~isempty(record) && numel(pt_ids)==size(record.case_4chs, 1))
        case_4ch = record.case_4chs{pt};
        case_2ch = record.case_2chs{pt};
        case_3ch = record.case_3chs{pt};
        case_sax = record.case_saxs{pt};
        case_lvot = record.case_lvots{pt};
        case_aov = record.case_aovs{pt};
        case_rv = record.case_rvs{pt};
        case_aorticarch = record.case_aorticarchs{pt};
    else
        case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '4ch');
        case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '2ch');
        case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '3ch');
        case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'sa');
        case_lvot = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'lvot');
        case_aov = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'aov');
        case_rv = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'rv');
        case_aorticarch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'aorticarch');
    end
    
    if(size(case_4ch, 1)+size(case_2ch, 1)+size(case_3ch, 1)+size(case_lvot, 1)+size(case_aov, 1)+size(case_rv, 1)==0)
        continue;
    end
    
    try
        study_date = case_4ch.study_dates(end,:);
    catch
        try
            study_date = case_2ch.study_dates(end,:);
        catch
            try
                study_date = case_sax.study_dates(end,:);
            catch
                try
                    study_date = case_3ch.study_dates(end,:);
                catch
                    try
                        study_date = case_lvot.study_dates(end,:);
                    catch
                        try
                            study_date = case_aov.study_dates(end,:);
                        catch
                            study_date = case_rv.study_dates(end,:);
                        end
                    end
                end
            end
        end
    end
           
    case_prefix = [study_date '_' pt_id];
    dst_dir = fullfile(aiDir, study_date, pt_id);
%     if(exist(dst_dir))
%         rmdir(dst_dir, 's');
%     end
%     
%     continue; 
%     

    num_cases = size(case_4ch, 1) + size(case_3ch, 1) + size(case_2ch, 1) + size(case_lvot, 1) + size(case_aov, 1) + size(case_rv, 1);
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id ' - ' study_date ' - ' num2str(num_cases)]);

    if(num_cases==0 && exist(dst_dir))
        rmdir(dst_dir, 's');
        continue;
    end
    if(~exist(dst_dir))
        mkdir(dst_dir);
    end
          
        % -----------------------------
        
        if(size(case_4ch, 1)>0) 
            process_one_view_all(dataDir, case_4ch, 'ch4', resDir, series_num, gmap_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag);
        end
        
        % -----------------------------
        
        if(size(case_3ch, 1)>0)
            process_one_view_all(dataDir, case_3ch, 'ch3', resDir, series_num, gmap_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag);
        end
        % -----------------------------
        
        if(size(case_2ch, 1)>0)
            process_one_view_all(dataDir, case_2ch, 'ch2', resDir, series_num, gmap_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag);
        end
        % -----------------------------
        
        if(size(case_lvot, 1)>0)
            process_one_view_all(dataDir, case_lvot, 'lvot', resDir, series_num, gmap_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag);
        end
        % -----------------------------
        
        if(size(case_aov, 1)>0)
            process_one_view_all(dataDir, case_aov, 'aov', resDir, series_num, gmap_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag);
        end
        % -----------------------------
        
        if(size(case_rv, 1)>0)
            process_one_view_all(dataDir, case_rv, 'rv', resDir, series_num, gmap_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag);
        end
        
        if(size(case_aorticarch, 1)>0)
            process_one_view_all(dataDir, case_rv, 'aorticarch', resDir, series_num, gmap_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag);
        end
        
        if(size(case_sax, 1)>0)
            process_one_sax_view_all(dataDir, case_sax, resDir, series_num, gmap_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag);
        end
    
    closeall
    closeall    
end
end

function process_one_view_all(dataDir, case_4ch, view_str, resDir, series_num, gmap_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag)
    for ii=1:size(case_4ch, 1)
        try
            fname = case_4ch.file_names{ii};
            case_4ch_dir = fullfile(resDir, case_4ch.study_dates(ii,:), case_4ch.file_names{ii});
            [path, sname, ext] = fileparts(case_4ch_dir); 
            sname= sname(~isspace(sname));
            case_prefix = [view_str '_' sname];

            dst_dir_4ch = fullfile(dst_dir, case_prefix);

            [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(fname);
            h5_name = fullfile(dataDir, study_dates, [fname '.h5']);

    %         if(check_processed && exist(fullfile(dst_dir, case_prefix, 'data.npy')) && exist(fullfile(dst_dir_4ch, 'gmap_slc_1.npy')) && exist(fullfile(dst_dir_4ch, 'images_for_gmap_imag.npy')))
    %             continue;
    %         end

            disp([case_prefix]);

            if(~check_processed || ~exist(fullfile(dst_dir, case_prefix, 'data.npy')))
                try
                    [gt_4ch, gt_h_4ch, acq_time, physio_time, attribs] = readGTPlusExportImageSeries_Squeeze(case_4ch_dir, series_num, 1);
                    gt_4ch = squeeze(gt_4ch);
                    gt_4ch = abs(gt_4ch);
                catch
                    disp(['Loading failed for ' case_4ch_dir ' - ' num2str(series_num)]);
                    continue
                end        

                process_one_view(view_str, pic_dir, case_prefix, ED, ES, dst_dir, gt_4ch, gt_h_4ch, dst_pixel_spacing, acq_time, physio_time, attribs, plot_flag);
            end

            if(~check_processed || ~exist(fullfile(dst_dir, case_prefix, 'gmap_acq.npy')) && gmap_num>=0)
                try
                    [gmap, gmap_h, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_4ch_dir, gmap_num, 1);
                    gmap = squeeze(gmap);
                    gmap = abs(gmap);

                    process_one_view_gmap(pic_dir, case_prefix, dst_dir, gmap, gmap_h, dst_pixel_spacing, plot_flag) 
                catch
                    disp(['Loading failed for ' case_4ch_dir ' - ' num2str(gmap_num)]);
                end      
            end

            if(~check_processed || ~exist(fullfile(dst_dir_4ch, 'ismrmrd_header.mat')))
                dset = ismrmrd.Dataset(h5_name, 'dataset');
                hdr = ismrmrd.xml.deserialize(dset.readxml);
                dset.close();

                save(fullfile(dst_dir_4ch, 'ismrmrd_header'), 'hdr');
            end
                
            debug_dir = fullfile(case_4ch_dir, 'DebugOutput');
            if(exist(debug_dir) & gmap_num>=0)                
                gmap_slices = [];
                accelFactor = [2 3 4 5];

                if(~check_processed || ~exist(fullfile(dst_dir_4ch, 'images_for_gmap_imag.npy')))
                    ori_images = readGTPlusExportData(fullfile(debug_dir, 'complexIm_afterRecon_afterPostProcessing'));
                    ori_images = squeeze(ori_images);
                    disp(['--> save into : ' dst_dir_4ch ]);
                    writeNPY(single(real(ori_images)), fullfile(dst_dir_4ch, 'images_for_gmap_real.npy'));
                    writeNPY(single(imag(ori_images)), fullfile(dst_dir_4ch, 'images_for_gmap_imag.npy'));            
                else
                    ori_images = complex(readNPY(fullfile(dst_dir_4ch, 'images_for_gmap_real.npy')), readNPY(fullfile(dst_dir_4ch, 'images_for_gmap_imag.npy')));
                end

                SLC = hdr.encoding(1).encodingLimits.slice.maximum + 1;
                for slc=1:SLC
                    gmap_file = fullfile(dst_dir_4ch, ['gmap_slc_' num2str(slc) '.npy']);
                    if(check_processed && exist(gmap_file))
                        continue;
                    end

                    acs_src = readGTPlusExportData(fullfile(debug_dir, ['acsSrc_n_0s_' num2str(slc-1)]));
                    acs_dst = readGTPlusExportData(fullfile(debug_dir, ['acsDst_n_0s_' num2str(slc-1)]));

                    coil_map = readGTPlusExportData(fullfile(debug_dir, ['coilMap_n_0s_' num2str(slc-1)]));

                    acs_src = squeeze(acs_src);
                    acs_dst = squeeze(acs_dst);
                    coil_map = squeeze(coil_map);

                    RO = size(coil_map, 1);
                    E1 = size(coil_map, 2);
                    srcCHA = size(acs_src, 3);
                    dstCHA = size(acs_dst, 3);

                    clear gmap
                    for af=1:numel(accelFactor)
                        kRO = 5;
                        kE1 = 4;
                        fitItself = 1;
                        thres = 5e-4;
                        [ker, convKer] = Matlab_gt_grappa_2d_calibrate(double(acs_src), double(acs_dst), kRO, kE1, accelFactor(af), fitItself, thres);

                        kIm = Matlab_gt_grappa_2d_compute_image_domain_kernel(double(convKer), RO, E1);

                    %     figure; imagescn(abs(kIm), [], [], [], 3)

                        [unmixing, gmap(:,:,af)] = Matlab_gt_grappa_2d_compute_unmxing_coeff(kIm, coil_map, accelFactor(af));
                    end

                    gmap = Matlab_gt_resize_2D_image(gmap, size(ori_images,1), size(ori_images,2), 5);

                    writeNPY(single(gmap), gmap_file);
                end            
            end
        catch
            if(exist(dst_dir_4ch))
                try
                    rmdir(dst_dir_4ch);
                catch
                end
            end
        end
    end
end

function process_one_view(view_str, pic_dir, case_prefix, ED, ES, dst_dir, gt_4ch, gt_h_4ch, dst_pixel_spacing, acq_time, physio_time, attribs, plot_flag)

    useMask = 1;
    boxFilterSize = 25;
    noisebackground = 30;
    thresRatioForNoise = 10;
    
    dst_dir_4ch = fullfile(dst_dir, case_prefix);
    mkdir(dst_dir_4ch);

    if(numel(size(gt_4ch))==4)
        gt_4ch = squeeze(gt_4ch(:,:, round(size(gt_4ch,3)/2),:));
    end

    [a_eigen,V,D] = KL_Eigenimage(gt_4ch);
%     if(isunix())
%         [scc, mask] = libMatlab_gt_surface_coil_correction(abs(double(a_eigen(:,:,end))), abs(double(a_eigen(:,:,end))), [], useMask, boxFilterSize, noisebackground, thresRatioForNoise);
%     else
        [scc, mask] = Matlab_gt_surface_coil_correction(abs(double(a_eigen(:,:,end))), abs(double(a_eigen(:,:,end))), [], useMask, boxFilterSize, noisebackground, thresRatioForNoise);
%     end
    
    ind = find(scc==inf);
    scc(ind) = 1;
    b = gt_4ch ./ repmat(scc, [1 1 size(gt_4ch,3)]);
    b(ind) = 0;
    gt_4ch_norm = b/max(b(:)) * max(gt_4ch(:));

    if(plot_flag)
        figure; imagescn(permute(gt_4ch_norm, [2 1 3]), [], [], [8], 3);
    end 

    RO = size(gt_4ch,1);
    E1 = size(gt_4ch,2);
    ps = max(gt_h_4ch(1,1,1).FOV)/max(RO,E1);

    new_RO = round(RO*ps/dst_pixel_spacing(1));
    new_E1 = round(E1*ps/dst_pixel_spacing(2));

    data = Matlab_gt_resize_2D_image(double(gt_4ch), new_RO, new_E1, 5);
    data_dst = Matlab_gt_resize_2D_image(double(gt_4ch_norm), new_RO, new_E1, 5);

    if(plot_flag)
        figure; imagescn(permute(data_dst, [2 1 3]), [], [], [8], 3);
    end

    gt_4ch = permute(gt_4ch, [2 1 3]);
    gt_4ch_norm = permute(gt_4ch_norm, [2 1 3]);

    data = permute(data, [2 1 3]);
    data_dst = permute(data_dst, [2 1 3]);
    size(data_dst)

    RO = size(data_dst,1);
    E1 = size(data_dst,2);
    PHS = size(data_dst,3);

%     data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), PHS);
%     data_dst = zpad(data_dst,max(RO, NN_RO), max(E1, NN_E1), PHS);

%     RO = size(data_dst,1);
%     E1 = size(data_dst,2);
%     PHS = size(data_dst,3);

    if(PHS<ES)
%         rmdir(dst_dir, 's');
        return;
    end

%     [mask, ro_s, ro_e, e1_s, e1_e, h_mask] = detect_heart_region_cine_roi(data_dst, 0.9, NN_RO, NN_E1, use_moco, regularization_hilbert_strength);
%     set(h_mask, 'visible', 'off');
% 
%     Cine_resized_training = data(ro_s:ro_e, e1_s:e1_e, :, :);
%     Cine_resized_training_norm = data_dst(ro_s:ro_e, e1_s:e1_e, :, :);
%     h2 = figure('visible', visible_status); imagescn(Cine_resized_training, [], [], [8], 3);
%     set(h2, 'visible', 'off')
%     h3 = figure('visible', visible_status); imagescn(Cine_resized_training_norm, [], [], [8], 3);
%     set(h3, 'visible', 'off')

%     header = CreateGtImageHeader(Cine_resized_training_norm);
%     Matlab_gt_write_analyze(single(Cine_resized_training), header, fullfile(dst_dir_4ch, 'Cine_resized_training'));
%     Matlab_gt_write_analyze(single(Cine_resized_training_norm), header, fullfile(dst_dir_4ch, 'Cine_resized_training_norm'));

    writeNPY(single(acq_time), fullfile(dst_dir_4ch, 'acq_time.npy'));
    writeNPY(single(physio_time), fullfile(dst_dir_4ch, 'physio_time.npy'));

    writeNPY(single(gt_4ch), fullfile(dst_dir_4ch, 'data_acq.npy'));
    writeNPY(single(gt_4ch_norm), fullfile(dst_dir_4ch, 'data_acq_norm.npy'));

%     writeNPY(single(Cine_resized_training), fullfile(dst_dir_4ch, 'Cine_resized_training.npy'));
%     writeNPY(single(Cine_resized_training_norm), fullfile(dst_dir_4ch, 'Cine_resized_training_norm.npy'));

    writeNPY(single(data), fullfile(dst_dir_4ch, 'data.npy'));
    writeNPY(single(data_dst), fullfile(dst_dir_4ch, 'data_norm.npy'));
%     writeNPY(single([ro_s, ro_e, e1_s, e1_e]), fullfile(dst_dir_4ch, 'roi.npy'));

%     saveas(h_mask, fullfile(dst_dir_4ch, 'roi'), 'jpg');
%     saveas(h2, fullfile(dst_dir_4ch, 'Cine_resized_training'), 'jpg');
%     saveas(h3, fullfile(dst_dir_4ch, 'Cine_resized_training_norm'), 'jpg');

    pixel_spacing = ps;
    gt_h = gt_h_4ch;
    save(fullfile(dst_dir_4ch, 'record'), 'gt_h', 'dst_pixel_spacing', 'pixel_spacing');
    save(fullfile(dst_dir_4ch, 'record_header'), 'gt_h', 'dst_pixel_spacing', 'pixel_spacing');

    ED_phs = data_dst(:,:,ED);
    ES_phs = data_dst(:,:,ES);

    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_' view_str '_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_' view_str '_' 'ES.jpg']), 'jpg');

    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(pic_dir, view_str, [case_prefix '_' view_str '_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(pic_dir, view_str, [case_prefix '_' view_str '_' 'ES.jpg']), 'jpg');

    writeNPY(single(ED_phs), fullfile(pic_dir, [view_str '_numpy'], [case_prefix '_' view_str '_' 'ED.npy']));
    writeNPY(single(ES_phs), fullfile(pic_dir, [view_str '_numpy'], [case_prefix '_' view_str '_' 'ES.npy']));

    ED_phs = data(:,:,ED);
    ES_phs = data(:,:,ES);

    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_ch4_original_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_ch4_original_' 'ES.jpg']), 'jpg');

    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(pic_dir, [view_str '_original'], [case_prefix '_' view_str '_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(pic_dir, [view_str '_original'], [case_prefix '_' view_str '_' 'ES.jpg']), 'jpg');

    writeNPY(single(ED_phs), fullfile(pic_dir, [view_str '_original_numpy'], [case_prefix '_' view_str '_' 'ED.npy']));
    writeNPY(single(ES_phs), fullfile(pic_dir, [view_str '_original_numpy'], [case_prefix '_' view_str '_' 'ES.npy']));
    
    save(fullfile(dst_dir_4ch, 'attributes'), 'attribs');
end

function process_one_sax_view_all(dataDir, case_sax, resDir, series_num, gmap_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag)
    for ii=1:size(case_sax, 1)        
        try
            fname = case_sax.file_names{ii};
            case_sax_dir = fullfile(resDir, case_sax.study_dates(ii,:), case_sax.file_names{ii});
            [path, sname, ext] = fileparts(case_sax.file_names{ii}); 
            sname= sname(~isspace(sname));
            case_prefix = ['sax_' sname];

            [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(fname);
            h5_name = fullfile(dataDir, study_dates, [fname '.h5']);

            dst_dir_sax = fullfile(dst_dir, case_prefix);

    %         if(check_processed && exist(fullfile(dst_dir, case_prefix, 'data.npy')) && exist(fullfile(dst_dir_sax, 'gmap_slc_1.npy')) && exist(fullfile(dst_dir_4ch, 'images_for_gmap_imag.npy')))
    %             continue;
    %         end

            SLC = case_sax.headers{ii}.encoding.encodingLimits.slice.maximum +1;

            accelFactor = [2 3 4 5];

            case_sax_dir = fullfile(resDir, case_sax.study_dates(ii,:), case_sax.file_names{ii});

            if(~check_processed || ~exist(fullfile(dst_dir, case_prefix, 'data.npy')))
                multi_slice_mode = 1;
                try   
                    [gt_sax_all, gt_h_sax_all, acq_time, physio_time, attribs] = readGTPlusExportImageSeries_Squeeze(case_sax_dir, series_num, 1);
                    gt_sax_all = squeeze(gt_sax_all);
                    multi_slice_mode = 0;
                catch
                    disp(['Loading failed for ' case_sax_dir ' - ' num2str(series_num)]);
                    continue;
                end

                if(size(gt_sax_all,4)==1 && SLC>1)
                    multi_slice_mode = 1;
                end

                if(SLC==1)
                    gt_sax_all = reshape(gt_sax_all, [size(gt_sax_all,1), size(gt_sax_all,2), 1, size(gt_sax_all,3)]);
                end

                process_one_sax_view(pic_dir, case_prefix, ED, ES, dst_dir, gt_sax_all, gt_h_sax_all, dst_pixel_spacing, acq_time, physio_time, attribs, multi_slice_mode, plot_flag);
            end

            if(~check_processed || ~exist(fullfile(dst_dir, case_prefix, 'gmap_acq.npy')) & gmap_num>=0)
                try
                    [gmap, gmap_h, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_sax_dir, gmap_num, 1);
                    gmap = squeeze(gmap);
                    gmap = abs(gmap);

                    process_one_view_gmap(pic_dir, case_prefix, dst_dir, gmap, gmap_h, dst_pixel_spacing, plot_flag) 
                catch
                    disp(['Loading failed for ' case_sax_dir ' - ' num2str(gmap_num)]);
                end    
            end

            debug_dir = fullfile(case_sax_dir, 'DebugOutput');
            if(~exist(debug_dir))
                return;
            end
            
            dset = ismrmrd.Dataset(h5_name, 'dataset');
            hdr = ismrmrd.xml.deserialize(dset.readxml);
            dset.close();

            save(fullfile(dst_dir_sax, 'ismrmrd_header'), 'hdr');

            gmap_slices = [];

            if(~check_processed || ~exist(fullfile(dst_dir_sax, 'images_for_gmap_imag.npy')) & gmap_num>=0)
                ori_images = readGTPlusExportData(fullfile(debug_dir, 'complexIm_afterRecon_afterPostProcessing'));
                ori_images = squeeze(ori_images);
                disp(['--> save into : ' dst_dir_sax ]);
                writeNPY(single(real(ori_images)), fullfile(dst_dir_sax, 'images_for_gmap_real.npy'));
                writeNPY(single(imag(ori_images)), fullfile(dst_dir_sax, 'images_for_gmap_imag.npy'));
            else
                if(gmap_num>=0)
                    ori_images = complex(readNPY(fullfile(dst_dir_sax, 'images_for_gmap_real.npy')), readNPY(fullfile(dst_dir_sax, 'images_for_gmap_imag.npy')));
                end
            end

            if(gmap_num>=0)
                SLC = hdr.encoding(1).encodingLimits.slice.maximum + 1;
                for slc=1:SLC
                    gmap_file = fullfile(dst_dir_sax, ['gmap_slc_' num2str(slc) '.npy']);
                    if(check_processed && exist(gmap_file))
                        continue;
                    end

                    acs_src = readGTPlusExportData(fullfile(debug_dir, ['acsSrc_n_0s_' num2str(slc-1)]));
                    acs_dst = readGTPlusExportData(fullfile(debug_dir, ['acsDst_n_0s_' num2str(slc-1)]));

                    coil_map = readGTPlusExportData(fullfile(debug_dir, ['coilMap_n_0s_' num2str(slc-1)]));

                    acs_src = squeeze(acs_src);
                    acs_dst = squeeze(acs_dst);
                    coil_map = squeeze(coil_map);

                    RO = size(coil_map, 1);
                    E1 = size(coil_map, 2);
                    srcCHA = size(acs_src, 3);
                    dstCHA = size(acs_dst, 3);

                    clear gmap
                    for af=1:numel(accelFactor)
                        kRO = 5;
                        kE1 = 4;
                        fitItself = 1;
                        thres = 5e-4;
                        [ker, convKer] = Matlab_gt_grappa_2d_calibrate(double(acs_src), double(acs_dst), kRO, kE1, accelFactor(af), fitItself, thres);

                        kIm = Matlab_gt_grappa_2d_compute_image_domain_kernel(double(convKer), RO, E1);

                    %     figure; imagescn(abs(kIm), [], [], [], 3)

                        [unmixing, gmap(:,:,af)] = Matlab_gt_grappa_2d_compute_unmxing_coeff(kIm, coil_map, accelFactor(af));
                    end

                    gmap = Matlab_gt_resize_2D_image(gmap, size(ori_images,1), size(ori_images,2), 5);

                    disp(['--> save gmap  ' gmap_file]);
                    writeNPY(single(gmap), gmap_file);
                end 
            end
        catch
            if(exist(dst_dir_sax))
                try
                    rmdir(dst_dir_sax);
                catch
                end
            end
        end
    end
end

function process_one_sax_view(pic_dir, case_prefix, ED, ES, dst_dir, gt_sax_all, gt_h_sax_all, dst_pixel_spacing, acq_time, physio_time, attribs, multi_slice_mode, plot_flag) 
    dst_dir_sax = fullfile(dst_dir, case_prefix);
    mkdir(dst_dir_sax);

    useMask = 1;
    boxFilterSize = 25;
    noisebackground = 30;
    thresRatioForNoise = 10;
    
    try
        gt_sax_norm = gt_sax_all;
        if(multi_slice_mode)
            gt_sax_norm = [];
            gt_sax_all = [];
            clear gt_h_sax_all
        end

        SLC = size(gt_sax_all, 3);
        
        for slc=1:SLC

            disp(['process slice ' num2str(slc) ' ... ']);

            if(multi_slice_mode)
                [gt_sax, gt_h_sax, acq_time, physio_time, attribs] = readGTPlusExportImageSeries_Squeeze(case_sax_dir, (slc-1)*100, 1);
                gt_sax = squeeze(gt_sax);

                gt_h_sax = gt_h_sax(slc, :, :);

                if(size(gt_sax,4)>1)
                    gt_sax = squeeze(gt_sax(:,:,slc,:));
                end

                gt_h_sax_all(slc) = gt_h_sax(1);
            else
                gt_h_sax = gt_h_sax_all(slc, 1, 1);
                gt_sax = squeeze(gt_sax_all(:,:,slc,:));
            end

            [a_eigen,V,D] = KL_Eigenimage(gt_sax);
            [scc, mask] = Matlab_gt_surface_coil_correction(abs(double(a_eigen(:,:,end))), abs(double(a_eigen(:,:,end))), [], useMask, boxFilterSize, noisebackground, thresRatioForNoise);

            ind = find(scc==inf);
            scc(ind) = 1;
            b = gt_sax ./ repmat(scc, [1 1 size(gt_sax,3)]);
            b(ind) = 0;
            gt_sax_norm(:,:,slc,:) = b/max(b(:)) * max(gt_sax(:));
            gt_sax_all(:,:,slc,:) = gt_sax;
        end

        if(plot_flag)
            figure; imagescn(permute(gt_sax_norm, [2 1 3 4]), [], [], [8], 4);
        end
    catch
        return
    end   

    % resample sax
    RO = size(gt_sax_all,1)
    E1 = size(gt_sax_all,2)
    ps = max(gt_h_sax_all(1,1,1).FOV)/max(RO,E1)

    new_RO = round(RO*ps/dst_pixel_spacing(1));
    new_E1 = round(E1*ps/dst_pixel_spacing(2));

    if(mod(new_RO, 2)==1)
        new_RO = new_RO + 1;
    end
    if(mod(new_E1, 2)==1)
        new_E1 = new_E1 + 1;
    end
    
    data = Matlab_gt_resize_2D_image(double(gt_sax_all), new_RO, new_E1, 5);
    data_dst = Matlab_gt_resize_2D_image(double(gt_sax_norm), new_RO, new_E1, 5);

    if(plot_flag)
        figure; imagescn(permute(data_dst, [2 1 3 4]), [], [], [8], 4);
    end

    gt_sax_all = permute(gt_sax_all, [2 1 4 3]);
    gt_sax_norm = permute(gt_sax_norm, [2 1 4 3]);

    data = permute(data, [2 1 4 3]);
    size(data)
    data_dst = permute(data_dst, [2 1 4 3]);
    size(data_dst)

    RO = size(data_dst,1);
    E1 = size(data_dst,2);
    PHS = size(data_dst,3);
    SLC = size(data_dst, 4);

%         [mask, ro_s, ro_e, e1_s, e1_e, h_mask] = detect_heart_region_cine_roi(data_dst, 0.95, NN_RO, NN_E1, use_moco, regularization_hilbert_strength);
% 
%         Cine_resized_training = data(ro_s:ro_e, e1_s:e1_e, :, :);
%         Cine_resized_training_norm = data_dst(ro_s:ro_e, e1_s:e1_e, :, :);
%         h2 = figure('visible', visible_status); imagescn(Cine_resized_training, [], [], [8], 3);
%         set(h2, 'visible', 'off')
%         h3 = figure('visible', visible_status); imagescn(Cine_resized_training_norm, [], [], [8], 3);
%         set(h3, 'visible', 'off')

%         header = CreateGtImageHeader(Cine_resized_training_norm);
%         Matlab_gt_write_analyze(single(Cine_resized_training), header, fullfile(dst_dir_sax, 'Cine_resized_training'));
%         Matlab_gt_write_analyze(single(Cine_resized_training_norm), header, fullfile(dst_dir_sax, 'Cine_resized_training_norm'));
% 
%         writeNPY(single(Cine_resized_training), fullfile(dst_dir_sax, 'Cine_resized_training.npy'));
%         writeNPY(single(Cine_resized_training_norm), fullfile(dst_dir_sax, 'Cine_resized_training_norm.npy'));

    writeNPY(single(acq_time), fullfile(dst_dir_sax, 'acq_time.npy'));
    writeNPY(single(physio_time), fullfile(dst_dir_sax, 'physio_time.npy'));
    
    writeNPY(single(gt_sax_all), fullfile(dst_dir_sax, 'data_acq.npy'));
    writeNPY(single(gt_sax_norm), fullfile(dst_dir_sax, 'data_acq_norm.npy'));

    writeNPY(single(data), fullfile(dst_dir_sax, 'data.npy'));
    writeNPY(single(data_dst), fullfile(dst_dir_sax, 'data_norm.npy'));
%         writeNPY(single([ro_s, ro_e, e1_s, e1_e]), fullfile(dst_dir_sax, 'roi.npy'));

%         saveas(h_mask, fullfile(dst_dir_sax, 'roi'), 'jpg');
%         saveas(h2, fullfile(dst_dir_sax, 'Cine_resized_training'), 'jpg');
%         saveas(h3, fullfile(dst_dir_sax, 'Cine_resized_training_norm'), 'jpg');

    pixel_spacing = ps;
    save(fullfile(dst_dir_sax, 'record'), 'gt_h_sax', 'dst_pixel_spacing', 'pixel_spacing');
    gt_h_sax = squeeze(gt_h_sax_all);

    save(fullfile(dst_dir_sax, 'attributes'), 'attribs');
    
    ED_phs = squeeze(data_dst(:,:,ED,:));
    ES_phs = squeeze(data_dst(:,:,ES,:));

    for slc=1:SLC
        im = ED_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(dst_dir_sax, [case_prefix '_sax_' 'ED_' num2str(slc) '.jpg']), 'jpg');
        im = ES_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(dst_dir_sax, [case_prefix '_sax_' 'ES_' num2str(slc) '.jpg']), 'jpg');
    end

    for slc=1:SLC
        im = ED_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(pic_dir, 'sax', [case_prefix '_sax_' 'ED_' num2str(slc) '.jpg']), 'jpg');
        im = ES_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(pic_dir, 'sax', [case_prefix '_sax_' 'ES_' num2str(slc) '.jpg']), 'jpg');

        writeNPY(single(ED_phs(:,:,slc)), fullfile(pic_dir, 'sax_numpy', [case_prefix '_sax_' 'ED_' num2str(slc) '.npy']));
        writeNPY(single(ES_phs(:,:,slc)), fullfile(pic_dir, 'sax_numpy', [case_prefix '_sax_' 'ES_' num2str(slc) '.npy']));
    end

    ED_phs = squeeze(data(:,:,ED,:));
    ES_phs = squeeze(data(:,:,ES,:));

    for slc=1:SLC
        im = ED_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(dst_dir_sax, [case_prefix '_sax_original_' 'ED_' num2str(slc) '.jpg']), 'jpg');
        im = ES_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(dst_dir_sax, [case_prefix '_sax_original_' 'ES_' num2str(slc) '.jpg']), 'jpg');
    end

    for slc=1:SLC
        im = ED_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(pic_dir, 'sax_original', [case_prefix '_sax_' 'ED_' num2str(slc) '.jpg']), 'jpg');
        im = ES_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(pic_dir, 'sax_original', [case_prefix '_sax_' 'ES_' num2str(slc) '.jpg']), 'jpg');

        writeNPY(single(ED_phs(:,:,slc)), fullfile(pic_dir, 'sax_original_numpy', [case_prefix '_sax_' 'ED_' num2str(slc) '.npy']));
        writeNPY(single(ES_phs(:,:,slc)), fullfile(pic_dir, 'sax_original_numpy', [case_prefix '_sax_' 'ES_' num2str(slc) '.npy']));
    end
   
    closeall
    closeall
end

% function process_one_view(view_str, pic_dir, case_prefix, ED, ES, dst_dir, gt_4ch, gt_h_4ch, dst_pixel_spacing, plot_flag)
% 
%     useMask = 1;
%     boxFilterSize = 25;
%     noisebackground = 30;
%     thresRatioForNoise = 10;
%     
%     dst_dir_4ch = fullfile(dst_dir, case_prefix);
%     mkdir(dst_dir_4ch);
% 
%     if(numel(size(gt_4ch))==4)
%         gt_4ch = squeeze(gt_4ch(:,:, round(size(gt_4ch,3)/2),:));
%     end
% 
%     [a_eigen,V,D] = KL_Eigenimage(gt_4ch);
% %     if(isunix())
% %         [scc, mask] = libMatlab_gt_surface_coil_correction(abs(double(a_eigen(:,:,end))), abs(double(a_eigen(:,:,end))), [], useMask, boxFilterSize, noisebackground, thresRatioForNoise);
% %     else
%         [scc, mask] = Matlab_gt_surface_coil_correction(abs(double(a_eigen(:,:,end))), abs(double(a_eigen(:,:,end))), [], useMask, boxFilterSize, noisebackground, thresRatioForNoise);
% %     end
%     
%     ind = find(scc==inf);
%     scc(ind) = 1;
%     b = gt_4ch ./ repmat(scc, [1 1 size(gt_4ch,3)]);
%     b(ind) = 0;
%     gt_4ch_norm = b/max(b(:)) * max(gt_4ch(:));
% 
%     if(plot_flag)
%         figure; imagescn(permute(gt_4ch_norm, [2 1 3]), [], [], [8], 3);
%     end 
% 
%     RO = size(gt_4ch,1);
%     E1 = size(gt_4ch,2);
%     ps = max(gt_h_4ch(1,1,1).FOV)/max(RO,E1);
% 
%     new_RO = round(RO*ps/dst_pixel_spacing(1));
%     new_E1 = round(E1*ps/dst_pixel_spacing(2));
% 
%     data = Matlab_gt_resize_2D_image(double(gt_4ch), new_RO, new_E1, 5);
%     data_dst = Matlab_gt_resize_2D_image(double(gt_4ch_norm), new_RO, new_E1, 5);
% 
%     if(plot_flag)
%         figure; imagescn(permute(data_dst, [2 1 3]), [], [], [8], 3);
%     end
% 
%     gt_4ch = permute(gt_4ch, [2 1 3]);
%     gt_4ch_norm = permute(gt_4ch_norm, [2 1 3]);
% 
%     data = permute(data, [2 1 3]);
%     data_dst = permute(data_dst, [2 1 3]);
%     size(data_dst)
% 
%     RO = size(data_dst,1);
%     E1 = size(data_dst,2);
%     PHS = size(data_dst,3);
% 
% %     data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), PHS);
% %     data_dst = zpad(data_dst,max(RO, NN_RO), max(E1, NN_E1), PHS);
% 
% %     RO = size(data_dst,1);
% %     E1 = size(data_dst,2);
% %     PHS = size(data_dst,3);
% 
%     if(PHS<ES)
% %         rmdir(dst_dir, 's');
%         return;
%     end
% 
% %     [mask, ro_s, ro_e, e1_s, e1_e, h_mask] = detect_heart_region_cine_roi(data_dst, 0.9, NN_RO, NN_E1, use_moco, regularization_hilbert_strength);
% %     set(h_mask, 'visible', 'off');
% % 
% %     Cine_resized_training = data(ro_s:ro_e, e1_s:e1_e, :, :);
% %     Cine_resized_training_norm = data_dst(ro_s:ro_e, e1_s:e1_e, :, :);
% %     h2 = figure('visible', visible_status); imagescn(Cine_resized_training, [], [], [8], 3);
% %     set(h2, 'visible', 'off')
% %     h3 = figure('visible', visible_status); imagescn(Cine_resized_training_norm, [], [], [8], 3);
% %     set(h3, 'visible', 'off')
% 
% %     header = CreateGtImageHeader(Cine_resized_training_norm);
% %     Matlab_gt_write_analyze(single(Cine_resized_training), header, fullfile(dst_dir_4ch, 'Cine_resized_training'));
% %     Matlab_gt_write_analyze(single(Cine_resized_training_norm), header, fullfile(dst_dir_4ch, 'Cine_resized_training_norm'));
% 
%     writeNPY(single(gt_4ch), fullfile(dst_dir_4ch, 'data_acq.npy'));
%     writeNPY(single(gt_4ch_norm), fullfile(dst_dir_4ch, 'data_acq_norm.npy'));
% 
% %     writeNPY(single(Cine_resized_training), fullfile(dst_dir_4ch, 'Cine_resized_training.npy'));
% %     writeNPY(single(Cine_resized_training_norm), fullfile(dst_dir_4ch, 'Cine_resized_training_norm.npy'));
% 
%     writeNPY(single(data), fullfile(dst_dir_4ch, 'data.npy'));
%     writeNPY(single(data_dst), fullfile(dst_dir_4ch, 'data_norm.npy'));
% %     writeNPY(single([ro_s, ro_e, e1_s, e1_e]), fullfile(dst_dir_4ch, 'roi.npy'));
% 
% %     saveas(h_mask, fullfile(dst_dir_4ch, 'roi'), 'jpg');
% %     saveas(h2, fullfile(dst_dir_4ch, 'Cine_resized_training'), 'jpg');
% %     saveas(h3, fullfile(dst_dir_4ch, 'Cine_resized_training_norm'), 'jpg');
% 
%     pixel_spacing = ps;
%     gt_h = gt_h_4ch;
%     save(fullfile(dst_dir_4ch, 'record'), 'gt_h', 'dst_pixel_spacing', 'pixel_spacing');
%     save(fullfile(dst_dir_4ch, 'record_header'), 'gt_h', 'dst_pixel_spacing', 'pixel_spacing');
% 
%     ED_phs = data_dst(:,:,ED);
%     ES_phs = data_dst(:,:,ES);
% 
%     imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_' view_str '_' 'ED.jpg']), 'jpg');
%     imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_' view_str '_' 'ES.jpg']), 'jpg');
% 
%     imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(pic_dir, view_str, [case_prefix '_' view_str '_' 'ED.jpg']), 'jpg');
%     imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(pic_dir, view_str, [case_prefix '_' view_str '_' 'ES.jpg']), 'jpg');
% 
%     writeNPY(single(ED_phs), fullfile(pic_dir, [view_str '_numpy'], [case_prefix '_' view_str '_' 'ED.npy']));
%     writeNPY(single(ES_phs), fullfile(pic_dir, [view_str '_numpy'], [case_prefix '_' view_str '_' 'ES.npy']));
% 
%     ED_phs = data(:,:,ED);
%     ES_phs = data(:,:,ES);
% 
%     imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_ch4_original_' 'ED.jpg']), 'jpg');
%     imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_ch4_original_' 'ES.jpg']), 'jpg');
% 
%     imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(pic_dir, [view_str '_original'], [case_prefix '_' view_str '_' 'ED.jpg']), 'jpg');
%     imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(pic_dir, [view_str '_original'], [case_prefix '_' view_str '_' 'ES.jpg']), 'jpg');
% 
%     writeNPY(single(ED_phs), fullfile(pic_dir, [view_str '_original_numpy'], [case_prefix '_' view_str '_' 'ED.npy']));
%     writeNPY(single(ES_phs), fullfile(pic_dir, [view_str '_original_numpy'], [case_prefix '_' view_str '_' 'ES.npy']));
% end

% function process_one_sax_view_all(case_sax, resDir, series_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag)
%     for ii=1:size(case_sax, 1)
%         
%         case_sax_dir = fullfile(resDir, case_sax.study_dates(ii,:), case_sax.file_names{ii});
%         [path, sname, ext] = fileparts(case_sax.file_names{ii}); 
%         sname= sname(~isspace(sname));
%         case_prefix = ['sax_' sname];
% 
%         if(check_processed && exist(fullfile(dst_dir, case_prefix, 'data.npy')))
%             continue;
%         end
%         
%         SLC = case_sax.headers{ii}.encoding.encodingLimits.slice.maximum +1;
%         
%         case_sax_dir = fullfile(resDir, case_sax.study_dates(ii,:), case_sax.file_names{ii})       
%         multi_slice_mode = 1;
%         try   
%             [gt_sax_all, gt_h_sax_all, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_sax_dir, series_num, 1);
%             gt_sax_all = squeeze(gt_sax_all);
%             multi_slice_mode = 0;
%         catch
%             disp(['Loading failed for ' case_sax_dir ' - ' num2str(series_num)]);
%             continue;
%         end
% 
%         if(size(gt_sax_all,4)==1 && SLC>1)
%             multi_slice_mode = 1;
%         end
%             
%         disp([case_prefix]);
%         
%         if(SLC==1)
%             gt_sax_all = reshape(gt_sax_all, [size(gt_sax_all,1), size(gt_sax_all,2), 1, size(gt_sax_all,3)]);
%         end
%         
%         process_one_sax_view(pic_dir, case_prefix, ED, ES, dst_dir, gt_sax_all, gt_h_sax_all, dst_pixel_spacing, multi_slice_mode, plot_flag);
%     end
% end

function process_one_view_gmap(pic_dir, case_prefix, dst_dir, gt_sax_all, gt_h_sax_all, dst_pixel_spacing, plot_flag) 
    dst_dir_sax = fullfile(dst_dir, case_prefix);
    mkdir(dst_dir_sax);
   
    SLC = size(gt_sax_all, 3);
  
    % resample sax
    RO = size(gt_sax_all,1)
    E1 = size(gt_sax_all,2)
    ps = max(gt_h_sax_all(1,1,1).FOV)/max(RO,E1)

    new_RO = round(RO*ps/dst_pixel_spacing(1));
    new_E1 = round(E1*ps/dst_pixel_spacing(2));

    if(mod(new_RO, 2)==1)
        new_RO = new_RO + 1;
    end
    if(mod(new_E1, 2)==1)
        new_E1 = new_E1 + 1;
    end
    
    data = Matlab_gt_resize_2D_image(double(gt_sax_all), new_RO, new_E1, 5);
    data = permute(data, [2 1 4 3]);
    gt_sax_all = permute(gt_sax_all, [2 1 4 3]);
    size(data)

    writeNPY(single(data), fullfile(dst_dir_sax, 'gmap.npy'));
    writeNPY(single(gt_sax_all), fullfile(dst_dir_sax, 'gmap_acq.npy'));
end