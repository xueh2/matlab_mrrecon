function apply_ai_cmr_view_detection_T1_MOLLI(resDir, aiDir, pt_ids, files_record_picked, cmr_view_model, script_name, checkprocessed)
% apply_ai_cmr_view_detection_T1_MOLLI(resDir, aiDir, pt_ids, files_record_picked, cmr_view_model, script_name, checkprocessed)

if(isunix())
    ext = '.sh';
else
    ext = '.bat';
end

% lax_scripts
if(isunix())
    command_pt = ['python3 '];
else
    command_pt = ['python '];
end
            
fid_cmr_view = fopen(fullfile(aiDir, [script_name '_cmr_view_detection' ext]), 'a+');

set_env_apply_ai(fid_cmr_view);

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    
    case_pre = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'pre');
    case_post = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'post');

    if(numel(case_pre)==0 || numel(case_post)==0)
        continue;
    end
    
    if(case_pre.headers{1}.subjectInformation.patientGender == 'O')
        continue;
    end        
           
    dst_dir = fullfile(aiDir, case_pre.study_dates(end,:), pt_id);
       
    for tt=1:2
        case_used = case_pre;
        if(tt==2)
            case_used = case_post;
        end

        for k=1:size(case_used, 1)

            t1_map_dir = fullfile(dst_dir, [case_used.file_names{k} '__t1_map']);
           
            if(exist(fullfile(t1_map_dir, 'data.npy')))
                contourDir = fullfile(dst_dir, [case_used.file_names{k} '__res_ai']);
                mkdir(contourDir);
            
                view_file = fullfile(contourDir, ['CMR_view_res_t1_map.npy']);
                probs_file = fullfile(contourDir, ['CMR_view_probs_t1_map.npy']);
                command = [command_pt ' cmr_view_detection.py --input ' fullfile(t1_map_dir, 'data.npy') ... 
                    ' --output ' view_file ' --prob ' probs_file ...
                    ' --model ' cmr_view_model];

                if(~checkprocessed || ~exist(view_file))
                    command
                    fprintf(fid_cmr_view, '%s\n', command);
                end
            end            
        end
    end
end

fclose(fid_cmr_view);
