function prepare_ai_perfusion_training_cmr_landmark(dataDir, resDir, aiDir, pt_ids, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_num, check_processed, plot_flag)
% prepare_ai_perfusion_training_cmr_landmark(dataDir, resDir, aiDir, pt_ids, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_num, check_processed, plot_flag)

mkdir(aiDir)

pic_dir = fullfile(aiDir, 'jpg_pics');
mkdir(pic_dir)
mkdir(fullfile(pic_dir, 'stress'))
mkdir(fullfile(pic_dir, 'rest'))
mkdir(fullfile(pic_dir, 'stress_numpy'))
mkdir(fullfile(pic_dir, 'rest_numpy'))
mkdir(fullfile(pic_dir, 'stress_pd_numpy'))
mkdir(fullfile(pic_dir, 'rest_pd_numpy'))
mkdir(fullfile(pic_dir, 'stress_pd_image_numpy'))
mkdir(fullfile(pic_dir, 'rest_pd_image_numpy'))

mkdir(fullfile(pic_dir, 'stress_pd_image'))
mkdir(fullfile(pic_dir, 'rest_pd_image'))

mkdir(fullfile(pic_dir, 'stress_eigen_image'))
mkdir(fullfile(pic_dir, 'rest_eigen_image'))
mkdir(fullfile(pic_dir, 'stress_eigen_image_numpy'))
mkdir(fullfile(pic_dir, 'rest_eigen_image_numpy'))

mkdir(fullfile(pic_dir, 'stress_ori_numpy'))
mkdir(fullfile(pic_dir, 'rest_ori_numpy'))
mkdir(fullfile(pic_dir, 'stress_ori_3T_numpy'))
mkdir(fullfile(pic_dir, 'rest_ori_3T_numpy'))

visible_status = 'on';
if(~plot_flag)
    visible_status = 'off';
end

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    
    case_stress = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'stress');
    case_rest = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'rest');

    if(numel(case_stress)==0 || numel(case_rest)==0 )
        continue;
    end
        
    case_prefix = [case_stress.study_dates(end,:) '_' pt_id];
    dst_dir = fullfile(aiDir, case_stress.study_dates(end,:), pt_id);
    
    stress_processed = 0;
    rest_processed = 0;
   
    if(check_processed)                
        if(exist(fullfile(pic_dir, 'stress_eigen_image_numpy', [case_prefix '_eigen_slc_3.npy'])))
            stress_processed = 1;
        end
        
        if(exist(fullfile(pic_dir, 'rest_eigen_image_numpy', [case_prefix '_eigen_slc_3.npy'])))
            rest_processed = 1;
        end        
    end
    
    try
        SLC = case_stress.headers{end}.encoding(1).encodingLimits.slice.maximum +1;
        sax_ind = size(case_stress, 1);
        if(SLC~=3)
            continue
        end

        % -----------------------------
        
        if(~stress_processed)
            case_stress_dir = fullfile(resDir, case_stress.study_dates(end,:), case_stress.file_names{end})
            
            perf0 = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'perf_moco_for_CASignal_0'));
            perf1 = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'perf_moco_for_CASignal_1'));
            perf2 = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'perf_moco_for_CASignal_2'));
            
            if(exist(fullfile(case_stress_dir, 'DebugOutput', 'input_for_moco_0_MAG.hdr')))
                ori_perf0 = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'input_for_moco_0_MAG'));
                ori_perf1 = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'input_for_moco_1_MAG'));
                ori_perf2 = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'input_for_moco_2_MAG'));
            else
                ori_im = readGTPlusExportImageSeries_Squeeze(case_stress_dir, 103);
                ori_perf0 = squeeze(ori_im(:,:,1,:));
                ori_perf1 = squeeze(ori_im(:,:,2,:));
                ori_perf2 = squeeze(ori_im(:,:,3,:));
            end
            
            pd0 = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'PD_for_moco_row0'));
            pd1 = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'PD_for_moco_row1'));
            pd2 = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'PD_for_moco_row2'));
            
            AIF_AcqTimes_0 = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'AIF_AcqTimes_0'));
            dstAcqTimes_0 = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'dstAcqTimes_0'));
            aif_mask_AI_LV = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'aif_mask_AI_LV'));
            aif_mask_AI_RV = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'aif_mask_AI_RV'));
            aif_RV_mask_final = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'aif_mask_AI_RV'));
            aif_LV_mask_final = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'aif_mask_AI_LV'));
            aif_moco = analyze75read(fullfile(case_stress_dir, 'DebugOutput', 'aif_moco'));
            
            RO = size(aif_moco, 1);
            E1 = size(aif_moco, 2);
            N = size(aif_moco, 3);
            
            if(RO ~= size(aif_LV_mask_final, 1))
                aif_moco = permute(aif_moco, [2 1 3]);
            end
            
            RV_s = squeeze(mean(reshape(aif_moco .* repmat(aif_RV_mask_final, [1, 1, N]), [RO*E1, N]), 1));
            LV_s = squeeze(mean(reshape(aif_moco .* repmat(aif_LV_mask_final, [1, 1, N]), [RO*E1, N]), 1));
            
            [mRV, peak_RV_s] = max(RV_s);
            [mLV, peak_LV_s] = max(LV_s);
            
%             [gt_stress, gt_h_stress, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_stress_dir, series_num, 1);
%             gt_stress = squeeze(gt_stress);
%             gt_stress = gt_stress(:,:,:,7:end);
            
             gt_stress = permute(cat(4, perf0, perf1, perf2), [2 1 3 4]);
             gt_stress_ori = permute(cat(4, ori_perf0, ori_perf1, ori_perf2), [2 1 3 4]);
             gt_stress_pd = permute(cat(3, pd0(:,:,1), pd1(:,:,1), pd2(:,:,1)), [2 1 3 4]);
             
             fs = case_stress.headers{end}.acquisitionSystemInformation.systemFieldStrength_T;
             if(fs>2)
                 writeNPY(single(gt_stress_ori), fullfile(pic_dir, 'stress_ori_3T_numpy', [case_stress.file_names{end} '.npy']));
             else
                 writeNPY(single(gt_stress_ori), fullfile(pic_dir, 'stress_ori_numpy', [case_stress.file_names{end} '.npy']));
             end
        end
        
        % -----------------------------
        
        if(~rest_processed)
            case_rest_dir = fullfile(resDir, case_rest.study_dates(end,:), case_rest.file_names{end})       
%             [gt_rest, gt_h_rest, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_rest_dir, series_num, 1);
%             gt_rest = squeeze(gt_rest);
%             gt_rest = gt_rest(:,:,:,7:end);

            perf0 = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'perf_moco_for_CASignal_0'));
            perf1 = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'perf_moco_for_CASignal_1'));
            perf2 = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'perf_moco_for_CASignal_2'));
            
            if(exist(fullfile(case_rest_dir, 'DebugOutput', 'input_for_moco_0_MAG.hdr')))
                ori_perf0 = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'input_for_moco_0_MAG'));
                ori_perf1 = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'input_for_moco_1_MAG'));
                ori_perf2 = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'input_for_moco_2_MAG'));
            else
                ori_im = readGTPlusExportImageSeries_Squeeze(case_rest_dir, 103);
                ori_perf0 = squeeze(ori_im(:,:,1,:));
                ori_perf1 = squeeze(ori_im(:,:,2,:));
                ori_perf2 = squeeze(ori_im(:,:,3,:));
            end
            
            pd0 = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'PD_for_moco_row0'));
            pd1 = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'PD_for_moco_row1'));
            pd2 = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'PD_for_moco_row2'));
            
            AIF_AcqTimes_0 = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'AIF_AcqTimes_0'));
            dstAcqTimes_0 = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'dstAcqTimes_0'));
            aif_mask_AI_LV = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'aif_mask_AI_LV'));
            aif_mask_AI_RV = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'aif_mask_AI_RV'));
            aif_RV_mask_final = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'aif_RV_mask_final'));
            aif_LV_mask_final = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'aif_LV_mask_final'));
            aif_moco = analyze75read(fullfile(case_rest_dir, 'DebugOutput', 'aif_moco'));
            
            RO = size(aif_moco, 1);
            E1 = size(aif_moco, 2);
            N = size(aif_moco, 3);
            
            if(RO ~= size(aif_LV_mask_final, 1))
                aif_moco = permute(aif_moco, [2 1 3]);
            end
            
            RV_r = squeeze(mean(reshape(aif_moco .* repmat(aif_RV_mask_final, [1, 1, N]), [RO*E1, N]), 1));
            LV_r = squeeze(mean(reshape(aif_moco .* repmat(aif_LV_mask_final, [1, 1, N]), [RO*E1, N]), 1));
            
            [mRV, peak_RV_r] = max(RV_r);
            [mLV, peak_LV_r] = max(LV_r);   
            
            gt_rest = permute(cat(4, perf0, perf1, perf2), [2 1 3 4]);
            gt_rest_ori = permute(cat(4, ori_perf0, ori_perf1, ori_perf2), [2 1 3 4]);
            gt_rest_pd = permute(cat(3, pd0(:,:,1), pd1(:,:,1), pd2(:,:,1)), [2 1 3 4]);
            
            fs = case_rest.headers{end}.acquisitionSystemInformation.systemFieldStrength_T;
             if(fs>2)
                 writeNPY(single(gt_rest_ori), fullfile(pic_dir, 'rest_ori_3T_numpy', [case_rest.file_names{end} '.npy']));
             else
                 writeNPY(single(gt_rest_ori), fullfile(pic_dir, 'rest_ori_numpy', [case_rest.file_names{end} '.npy']));
             end
        end
        
        if(plot_flag)
            h1 = figure; imagescn(gt_stress, [], [1 SLC], [8], 3);
            h1 = figure; imagescn(gt_rest, [], [1 SLC], [8], 3);
            h1 = figure; imagescn(gt_stress_pd, [], [1 SLC], [8], 3);
            h1 = figure; imagescn(gt_rest_pd, [], [1 SLC], [8], 3);
        end
    catch
        continue;
    end
       
    try
        mkdir(dst_dir)

        % stress
        if(~stress_processed)
            make_perf_pics('stress', dst_dir, gt_stress, gt_stress_pd, case_stress, dst_pixel_spacing, plot_flag, RV_s, LV_s, case_prefix, pic_dir);
        end
        
        % rest
        if(~rest_processed)
            make_perf_pics('rest', dst_dir, gt_stress, gt_rest_pd, case_rest, dst_pixel_spacing, plot_flag, RV_r, LV_r, case_prefix, pic_dir);            
        end
    catch
        rmdir(dst_dir, 's');
        continue;
    end
        
    closeall
    closeall    
end

end

function make_perf_pics(data_role, dst_dir, gt_stress, gt_stress_pd, case_stress, dst_pixel_spacing, plot_flag, RV, LV, case_prefix, pic_dir)
    dst_dir_stress = fullfile(dst_dir, data_role);
    mkdir(dst_dir_stress);

    [mRV, peak_RV_s] = max(RV);
    [mLV, peak_LV_s] = max(LV);
            
    RO = size(gt_stress,1)
    E1 = size(gt_stress,2)
    N = size(gt_stress, 3);
    SLC = size(gt_stress, 4);
    ps = max(case_stress.headers{end}.encoding(1).reconSpace.fieldOfView_mm.x, case_stress.headers{end}.encoding(1).reconSpace.fieldOfView_mm.y)/max(RO,E1)

    %new_RO = round(RO*ps/dst_pixel_spacing(1));
    %new_E1 = round(E1*ps/dst_pixel_spacing(2));

    new_RO = 2*RO;
    new_E1 = 2*E1;
    
    data = Matlab_gt_resize_2D_image(double(gt_stress), new_RO, new_E1, 5);
    if(plot_flag)
        figure; imagescn(data, [], [], [8], 3);
    end

    data_pd = Matlab_gt_resize_2D_image(double(gt_stress_pd), new_RO, new_E1, 5);

    RO = size(data,1);
    E1 = size(data,2);
    REP = size(data,3);
    SLC = size(data,4);

    if(REP<24)
        rmdir(dst_dir, 's');
        return;
    end

    a_eigen = zeros(RO, E1, SLC);

    for slc=1:SLC
        %[ss,V,D] = KL_Eigenimage(squeeze(data(:,:,peak_RV_s-1:peak_LV_s, slc)));
        ss = std(data(:,:,:,slc), [], 3);
        a_eigen(:,:,slc) = ss(:,:,end);
    end
    a_eigen = abs(a_eigen);
    if(plot_flag)
        figure; imagescn(a_eigen, [], [], [8]);
    end

    header = CreateGtImageHeader(data);
    Matlab_gt_write_analyze(single(data), header, fullfile(dst_dir_stress, 'data'));

    writeNPY(single(gt_stress), fullfile(dst_dir_stress, 'data_acq.npy'));
    writeNPY(single(data), fullfile(dst_dir_stress, 'data.npy'));

    writeNPY(single(a_eigen), fullfile(dst_dir_stress, 'data_eigen.npy'));

    header = CreateGtImageHeader(data_pd);
    Matlab_gt_write_analyze(single(data_pd), header, fullfile(dst_dir_stress, 'pd'));
    writeNPY(single(gt_stress_pd), fullfile(dst_dir_stress, 'pd_acq.npy'));
    writeNPY(single(data_pd), fullfile(dst_dir_stress, 'pd.npy'));

    pixel_spacing = ps;
    protocol_header = case_stress.headers{end};
    save(fullfile(dst_dir_stress, 'record'), 'protocol_header', 'dst_pixel_spacing', 'pixel_spacing', 'RV', 'LV');
    save(fullfile(dst_dir_stress, 'record_header'), 'protocol_header', 'dst_pixel_spacing', 'pixel_spacing', 'RV', 'LV');
   
    for slc=1:SLC
        im = a_eigen(:,:,slc);
        I = uint8(255*im/max(im(:)));
        J = imadjust(I,stretchlim(I, [0.01, 0.99]),[]);
        imwrite(J, fullfile(dst_dir_stress, [case_prefix '_eigen_slc_' num2str(slc) '.jpg']), 'jpg');
        copyfile(fullfile(dst_dir_stress, [case_prefix '_eigen_slc_' num2str(slc) '.jpg']), fullfile(pic_dir, [data_role '_eigen_image']));
        writeNPY(single(im), fullfile(pic_dir, [data_role '_eigen_image_numpy'], [case_prefix '_eigen_slc_' num2str(slc) '.npy']));

        im = data_pd(:,:,slc);
        I = uint8(255*im/max(im(:)));
        J = imadjust(I,stretchlim(I),[]);
        imwrite(J, fullfile(dst_dir_stress, [case_prefix '_pd_slc_' num2str(slc) '.jpg']), 'jpg');
        copyfile(fullfile(dst_dir_stress, [case_prefix '_pd_slc_' num2str(slc) '.jpg']), fullfile(pic_dir, [data_role '_pd_image']));
        writeNPY(single(im), fullfile(pic_dir, [data_role '_pd_image_numpy'], [case_prefix '_pd_slc_' num2str(slc) '.npy']));
    end

    writeNPY(single(data), fullfile(pic_dir, [data_role '_numpy'], [case_prefix '.npy']));
    writeNPY(single(data_pd), fullfile(pic_dir, [data_role '_pd_numpy'], [case_prefix '.npy']));
end
