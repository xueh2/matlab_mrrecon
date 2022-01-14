function prepare_ai_cine_training_retro_AI_seg(dataDir, resDir, aiDir, pt_ids, files_record_picked, record)
% prepare_ai_cine_training_retro_AI_seg(dataDir, resDir, aiDir, pt_ids, files_record_picked, record)

mkdir(aiDir)

checkProcessed = 0;
sendDicom = 1;
startRemoteGT = 0;
configNames = {'GTPrep_2DT_RetroCine_GLS_Seg_AI_istore.xml'}
copy_debug_output = 1
copy_dicom_output = 0

gt_host = 'localhost'
GT_PORT = '9014'

if(isunix())
    pre_set_debug_folder = '/home/xueh2/Debug/DebugOutput_RetroCine'
    mkdir(pre_set_debug_folder)
else
    pre_set_debug_folder = []
end

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
    
    if(size(case_4ch, 1)==0)
        continue;
    end
    if(size(case_2ch, 1)==0)
        continue;
    end
    if(size(case_sax, 1)==0)
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

    num_cases = size(case_4ch, 1) + size(case_3ch, 1) + size(case_2ch, 1) + size(case_lvot, 1) + size(case_aov, 1) + size(case_rv, 1);
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id ' - ' study_date ' - ' num2str(num_cases)]);

    dst_dir = fullfile(aiDir, study_date, pt_id);
    if(num_cases==0 && exist(dst_dir))
        rmdir(dst_dir, 's');
        continue;
    end    
          
    % -----------------------------

    if(size(case_4ch, 1)>0) 
        
        case_res_dir = fullfile(dst_dir, case_4ch.file_names{end, 1});
        if(~exist(case_res_dir))
            files_all = {case_4ch.file_names{end, 1}};
            [tUsed, ignored, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, files_all, gt_host, resDir, checkProcessed, sendDicom, startRemoteGT, configNames, [], GT_PORT, copy_debug_output, copy_dicom_output, pre_set_debug_folder);
        end
    end

    % -----------------------------

%     if(size(case_3ch, 1)>0)
%         process_one_view_all(dataDir, case_3ch, 'ch3', resDir, series_num, gmap_num, pic_dir, ED, ES, dst_dir, dst_pixel_spacing, check_processed, plot_flag);
%     end
    % -----------------------------

    if(size(case_2ch, 1)>0)
        
        case_res_dir = fullfile(dst_dir, case_2ch.file_names{end, 1});
        if(~exist(case_res_dir))
            files_all = {case_2ch.file_names{end, 1}};
            [tUsed, ignored, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, files_all, gt_host, resDir, checkProcessed, sendDicom, startRemoteGT, configNames, [], GT_PORT, copy_debug_output, copy_dicom_output, pre_set_debug_folder);
        end
    end
    % -----------------------------
    

    if(size(case_sax, 1)>0)

        case_res_dir = fullfile(dst_dir, case_sax.file_names{end, 1});
        if(~exist(case_res_dir))
            files_all = {case_sax.file_names{end, 1}};
            [tUsed, ignored, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, files_all, gt_host, resDir, checkProcessed, sendDicom, startRemoteGT, configNames, [], GT_PORT, copy_debug_output, copy_dicom_output, pre_set_debug_folder);
        end
    end
    
%     save(fullfile(dst_dir, 'record.mat'), 'case_4ch', 'case_2ch', 'case_sax');
    
    debug_dst_dir = fullfile(dst_dir, 'ai_data');
    if(~exist(debug_dst_dir))
        mkdir(dst_dir);
        
        case_4ch_res_dir = fullfile(aiDir, study_date, case_4ch.file_names{end, 1});
        case_2ch_res_dir = fullfile(aiDir, study_date, case_2ch.file_names{end, 1});
        case_sax_res_dir = fullfile(aiDir, study_date, case_sax.file_names{end, 1});
        
        if(exist(case_4ch_res_dir))
            movefile(case_4ch_res_dir, dst_dir);
        end
        if(exist(case_2ch_res_dir))
            movefile(case_2ch_res_dir, dst_dir);
        end
        if(exist(case_sax_res_dir))
            movefile(case_sax_res_dir, dst_dir);
        end
        
        debug_dir = fullfile(dst_dir, case_sax.file_names{end, 1}, 'DebugOutput');
        
        try
            cd(debug_dir)
            Dicom_image_position = analyze75read('Dicom_image_position.hdr');
            Dicom_image_pixel_spacing = analyze75read('Dicom_image_pixel_spacing.hdr');
            Dicom_image_orientation = analyze75read('Dicom_image_orientation.hdr');
            I_sax_before_reporting = analyze75read('I_sax_before_reporting.hdr');
            I_ch2_before_reporting = analyze75read('I_ch2_before_reporting.hdr');
            I_ch4_before_reporting = analyze75read('I_ch4_before_reporting.hdr');
            I_sax_norm_before_reporting = analyze75read('I_sax_norm_before_reporting.hdr');
            I_ch2_norm_before_reporting = analyze75read('I_ch2_norm_before_reporting.hdr');
            I_ch4_norm_before_reporting = analyze75read('I_ch4_norm_before_reporting.hdr');
            sax_endo_mask_1mm = analyze75read('sax_endo_mask_1mm.hdr');
            sax_epi_mask_1mm = analyze75read('sax_epi_mask_1mm.hdr');
            ch2_pts_1mm = analyze75read('ch2_pts_1mm.hdr');
            ch4_pts_1mm = analyze75read('ch4_pts_1mm.hdr');
            metrics_1mm = analyze75read('metrics_1mm.hdr');
        catch
            continue
        end
                
        mkdir(debug_dst_dir);
        
        writeNPY(Dicom_image_position, fullfile(debug_dst_dir, 'Dicom_image_position.npy'));
        writeNPY(Dicom_image_pixel_spacing, fullfile(debug_dst_dir, 'Dicom_image_pixel_spacing.npy'));
        writeNPY(Dicom_image_orientation, fullfile(debug_dst_dir, 'Dicom_image_orientation.npy'));
        writeNPY(I_sax_before_reporting, fullfile(debug_dst_dir, 'I_sax_before_reporting.npy'));
        writeNPY(I_ch2_before_reporting, fullfile(debug_dst_dir, 'I_ch2_before_reporting.npy'));
        writeNPY(I_ch4_before_reporting, fullfile(debug_dst_dir, 'I_ch4_before_reporting.npy'));
        writeNPY(I_sax_norm_before_reporting, fullfile(debug_dst_dir, 'I_sax_norm_before_reporting.npy'));
        writeNPY(I_ch2_norm_before_reporting, fullfile(debug_dst_dir, 'I_ch2_norm_before_reporting.npy'));
        writeNPY(I_ch4_norm_before_reporting, fullfile(debug_dst_dir, 'I_ch4_norm_before_reporting.npy'));
        writeNPY(sax_endo_mask_1mm, fullfile(debug_dst_dir, 'sax_endo_mask_1mm.npy'));
        writeNPY(sax_epi_mask_1mm, fullfile(debug_dst_dir, 'sax_epi_mask_1mm.npy'));
        writeNPY(ch2_pts_1mm, fullfile(debug_dst_dir, 'ch2_pts_1mm.npy'));
        writeNPY(ch4_pts_1mm, fullfile(debug_dst_dir, 'ch4_pts_1mm.npy'));
        writeNPY(metrics_1mm, fullfile(debug_dst_dir, 'metrics_1mm.npy'));
    end
    cd(resDir)
    closeall
    closeall    
end
end
