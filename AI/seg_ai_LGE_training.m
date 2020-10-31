function seg_ai_LGE_training(resDir, aiDir, files_record_picked)
% seg_ai_LGE_training(resDir, aiDir, files_record_picked)

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
    if(~exist(dst_dir) | ~exist(fullfile(dst_dir, 'sax', 'Cine_resized_training.jpg')))
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
        case_2ch_dir = fullfile(resDir, case_2ch.study_dates(end,:), case_2ch.file_names{end})       
        case_sax_dir = fullfile(resDir, case_sax.study_dates(sax_ind,:), case_sax.file_names{sax_ind})       
    catch
        continue;
    end
    
    mkdir(dst_dir)

    % 4ch
    dst_dir_4ch = fullfile(dst_dir, 'ch4');
          
    % 2ch
    dst_dir_2ch = fullfile(dst_dir, 'ch2');
    
    %sax
    dst_dir_sax = fullfile(dst_dir, 'sax');
    data = readNPY(fullfile(dst_dir_sax, 'data.npy'));
    
    dst_dir_sax_seg = fullfile(dst_dir, 'sax_seg');
    mkdir(dst_dir_sax_seg);
    cd(dst_dir_sax_seg);
    
    figure;
    imagescn(data, [4096-300 4096+300], [], [12]);
    pause
    
    closeall
    closeall    
end
