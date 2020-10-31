function plot_ai_LGE_training(resDir, aiDir, files_record_picked, linewidth)
% plot_ai_LGE_training(resDir, aiDir, files_record_picked, linewidth)

pt_ids = unique(files_record_picked(:, 3));
load_endo = 1;
load_epi = 1;

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
    
    try
        SLC = case_sax.headers{end}.encoding.encodingLimits.slice.maximum +1;
        sax_ind = size(case_sax, 1);
        if(SLC<6)
            SLC = case_sax.headers{1}.encoding.encodingLimits.slice.maximum +1;            
            sax_ind = 1;
        end

        case_4ch_dir = fullfile(resDir, case_4ch.study_dates(end,:), case_4ch.file_names{end})       
        case_2ch_dir = fullfile(resDir, case_2ch.study_dates(end,:), case_2ch.file_names{end})       
        case_sax_dir = fullfile(resDir, case_sax.study_dates(end,:), case_sax.file_names{sax_ind})       
    catch
        continue;
    end
    
    % 4ch
    dst_dir_4ch = fullfile(dst_dir, 'ch4');
               
    % 2ch
    dst_dir_2ch = fullfile(dst_dir, 'ch2');
           
    %sax
    dst_dir_sax = fullfile(dst_dir, 'sax');
       
    contourDir = fullfile(dst_dir, 'sax_seg');
   
    if(exist(fullfile(contourDir, 'roi_new.mat')))    
        if(~exist(fullfile(contourDir, 'NN_contours.mat')))
            
            data = readNPY(fullfile(dst_dir_sax, 'data.npy'));
            roi = load(fullfile(contourDir, 'roi_new.mat'));
            plotFlag = 1;
            S = generate_seg_from_manual_roi(data, roi.ROI_info_table, plotFlag);
            
            endo = S.endo_mask;
            epi = S.epi_mask;
            
            save(fullfile(contourDir, 'NN_contours'), 'endo', 'epi', 'S');

            writeNPY(single(S.endo_mask), fullfile(contourDir, 'endo_mask.npy'));
            writeNPY(single(S.epi_mask), fullfile(contourDir, 'epi_mask.npy'));
            writeNPY(single(S.rvi_mask), fullfile(contourDir, 'rvi_mask.npy'));
        end
    end
end
