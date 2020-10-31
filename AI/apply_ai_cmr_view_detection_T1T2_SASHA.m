function apply_ai_cmr_view_detection_T1T2_SASHA(resDir, aiDir, pt_ids, files_record_picked, cmr_view_model, script_name, checkprocessed)
% apply_ai_cmr_view_detection_T1T2_SASHA(resDir, aiDir, pt_ids, files_record_picked, cmr_view_model, script_name, checkprocessed)

if(isunix())
    ext = '.sh';
else
    ext = '.bat';
end

fid_cmr_view = fopen(fullfile(aiDir, [script_name '_cmr_view_t2w' ext]), 'a+');
set_env_apply_ai(fid_cmr_view);

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    
    case_t1t2 = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'T1T2');

    if(numel(case_t1t2)==0)
        continue;
    end
    
    if(case_t1t2.headers{1}.subjectInformation.patientGender == 'O')
        continue;
    end
        
    dst_dir = fullfile(aiDir, case_t1t2.study_dates(end,:), pt_id);
       
    case_used = case_t1t2;
        
    for k=1:size(case_used, 1)
        contourDir = fullfile(dst_dir, [case_used.file_names{k} '__res_ai']);
        mkdir(contourDir);
        
        t2w_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2w_image']);
        if(~exist(t2w_image_dir))
            continue;
        end
            
        if(isunix())
            command_pt = ['python3 '];
        else
            command_pt = ['python '];
        end

        if(exist(fullfile(t2w_image_dir, 'data.npy')))
            view_file = fullfile(contourDir, 'CMR_view_res.npy');
            probs_file = fullfile(contourDir, 'CMR_view_probs.npy');
            command = [command_pt ' cmr_view_detection.py --input ' fullfile(t2w_image_dir, 'data.npy') ... 
                ' --output ' view_file ' --prob ' probs_file ...
                ' --model ' cmr_view_model];

            if(~checkprocessed || ~exist(view_file))
                command
                fprintf(fid_cmr_view, '%s\n', command);
            end
        end
    end
end

fclose(fid_cmr_view);
