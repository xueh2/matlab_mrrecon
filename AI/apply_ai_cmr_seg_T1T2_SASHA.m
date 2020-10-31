function apply_ai_cmr_seg_T1T2_SASHA(resDir, aiDir, pt_ids, files_record_picked, cmr_seg_model, script_name, checkprocessed)
% apply_ai_cmr_seg_T1T2_SASHA(resDir, aiDir, pt_ids, files_record_picked, cmr_seg_model, script_name, checkprocessed)

if(isunix())
    ext = '.sh';
else
    ext = '.bat';
end

fid = fopen(fullfile(aiDir, [script_name '_cmr_seg' ext]), 'a+');
set_env_apply_ai(fid);

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
        
        view_file = fullfile(contourDir, 'CMR_view_res.npy');
        if(~exist(view_file))
            continue;
        end
        view_res = readNPY(view_file);
        if(view_res(1)~=3)
            continue;
        end
        
        t1w_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t1w_image']);
        t2w_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2w_image']);
        pd_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__pd_image']);
        t1_map_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t1_map']);
        t2_map_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2_map']);
        if(~exist(t1w_image_dir))
            continue;
        end
        if(~exist(t2w_image_dir))
            continue;
        end
        if(~exist(pd_image_dir))
            continue;
        end
            
        if(isunix())
            command_pt = ['python3 '];
        else
            command_pt = ['python '];
        end

        if(exist(fullfile(t2w_image_dir, 'data.npy')))
            
            mask_file = fullfile(contourDir, 'T1T2_seg_mask.npy');
            probs_file = fullfile(contourDir, 'T1T2_seg_prob.npy');
            endo_file = fullfile(contourDir, 'T1T2_seg_endo.npy');
            epi_file = fullfile(contourDir, 'T1T2_seg_epi.npy');
            input_for_AI_file = fullfile(contourDir, 'T1T2_seg_input_for_AI.npy');
            
            command = [command_pt ' T1T2_seg_sax.py --t1 ' fullfile(t1_map_image_dir, 'data.npy') ... 
                ' --t2 ' fullfile(t2_map_image_dir, 'data.npy') ...
                ' --t1w ' fullfile(t1w_image_dir, 'data.npy') ...
                ' --t2w ' fullfile(t2w_image_dir, 'data.npy') ...
                ' --pd ' fullfile(pd_image_dir, 'data.npy') ...
                ' --output ' mask_file ' --prob ' probs_file ...
                ' --endo ' endo_file ' --epi ' epi_file ...
                ' --input_for_AI ' input_for_AI_file ' --model_POSE ' cmr_seg_model];

            if(~checkprocessed || ~exist(mask_file))
                command
                fprintf(fid, '%s\n', command);
            end
        end
    end
end

fclose(fid);
