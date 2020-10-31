function copy_training_data_ai_cine_training_retro_binning_rt(resDir, aiDir, trainDir_used, pt_ids, files_record_picked)
% copy_training_data_ai_cine_training_retro_binning_rt(resDir, aiDir, trainDir_used, pt_ids, files_record_picked)

load_endo = 1;
load_epi = 1;

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    
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
       
    contourDir = fullfile(dst_dir, 'sax_ai');
   
     if(exist(fullfile(dst_dir_sax, 'data_norm.npy')))    

%         Cine_resized_training = readNPY(fullfile(dst_dir_sax, 'Cine_resized_training_norm.npy'));
%         im_roi = readNPY(fullfile(dst_dir_sax, 'roi.npy'));
%         data = readNPY(fullfile(dst_dir_sax, 'data.npy'));
%         data_norm = readNPY(fullfile(dst_dir_sax, 'data_norm.npy'));
         
        train_data_dir = fullfile(trainDir_used, case_4ch.study_dates(end,:), pt_id);
        mkdir(train_data_dir);
                
        copyfile(fullfile(dst_dir_sax, 'data_norm.npy'), fullfile(train_data_dir, 'sax.npy'));
     end
end
