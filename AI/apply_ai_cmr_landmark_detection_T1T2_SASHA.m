function apply_ai_cmr_landmark_detection_T1T2_SASHA(resDir, aiDir, pt_ids, files_record_picked, sax_model, script_name, checkprocessed)
% apply_ai_cmr_landmark_detection_T1T2_SASHA(resDir, aiDir, pt_ids, files_record_picked, sax_model, script_name, checkprocessed)

if(isunix())
    ext = '.sh';
else
    ext = '.bat';
end

fid_sax = fopen(fullfile(aiDir, [script_name '_cmr_landmark_sax' ext]), 'a+');
fid_sax_t2_map = fopen(fullfile(aiDir, [script_name '_cmr_landmark_sax_t2_map' ext]), 'a+');
fid_sax_t2w = fopen(fullfile(aiDir, [script_name '_cmr_landmark_sax_t2w' ext]), 'a+');

set_env_apply_ai(fid_sax);
set_env_apply_ai(fid_sax_t2_map);
set_env_apply_ai(fid_sax_t2w);

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
        
    if(isunix())
        command = ['python3 '];
    else
        command = ['python '];
    end 
    
    for k=1:size(case_used, 1)
        contourDir = fullfile(dst_dir, [case_used.file_names{k} '__res_ai']);
        mkdir(contourDir);
        
        t2w_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2w_image']);
        
        if(exist(fullfile(t2w_image_dir, 'data.npy')))
            pts_file = fullfile(contourDir, ['SAX_AI_pts_t2w_image.npy']);
            probs_file = fullfile(contourDir, ['SAX_AI_probs_t2w_image.npy']);
            sax_command = [command 'cmr_landmark_detection.py --input ' fullfile(t2w_image_dir, 'data.npy') ... 
                ' --output ' pts_file ' --prob ' probs_file ...
                ' --model ' sax_model ' --batch_size 16' ' --lax 0 --use_3D 0 --smooth_pts 0 --RO 352 --E1 352 --cli_mode 1'];

            if(~checkprocessed || ~exist(pts_file))
                sax_command
                fprintf(fid_sax_t2w, '%s\n', sax_command);
            end
        end
            
        for slc=1:6
            t2_last_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2_last_image_slc' num2str(slc)]);
            t2_map_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2_map_slc' num2str(slc)]);

            if(~exist(t2_last_image_dir))
                if(slc==1)
                    t2_last_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2_last_image']);
                    t2_map_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2_map']);
                    
                    if(~exist(t2_last_image_dir))
                        continue;
                    end
                else
                    continue;
                end
            end
               
            % sax_scripts
            % -----------------------------
            
            if(exist(fullfile(t2_last_image_dir, 'data.npy')))
                pts_file = fullfile(contourDir, ['SAX_AI_pts_t2_last_image_slc' num2str(slc) '.npy']);
                probs_file = fullfile(contourDir, ['SAX_AI_probs_t2_last_image_slc' num2str(slc) '.npy']);
                sax_command = [command 'cmr_landmark_detection.py --input ' fullfile(t2_last_image_dir, 'data.npy') ... 
                    ' --output ' pts_file ' --prob ' probs_file ...
                    ' --model ' sax_model ' --batch_size 16' ' --lax 0 --use_3D 0 --smooth_pts 0 --RO 352 --E1 352 --cli_mode 1'];

                if(~checkprocessed || ~exist(pts_file))
                    sax_command
                    fprintf(fid_sax, '%s\n', sax_command);
                end
            end
            
            % -----------------------------
            
            if(exist(fullfile(t2_map_dir, 'data.npy')))
                pts_file = fullfile(contourDir, ['SAX_AI_pts_t2_map_slc' num2str(slc) '.npy']);
                probs_file = fullfile(contourDir, ['SAX_AI_probs_t2_map_slc' num2str(slc) '.npy']);
                sax_command = [command 'cmr_landmark_detection.py --input ' fullfile(t2_map_dir, 'data.npy') ... 
                    ' --output ' pts_file ' --prob ' probs_file ...
                    ' --model ' sax_model ' --batch_size 16' ' --lax 0 --use_3D 0 --smooth_pts 0 --RO 352 --E1 352 --cli_mode 1'];

                if(~checkprocessed || ~exist(pts_file))
                    sax_command
                    fprintf(fid_sax_t2_map, '%s\n', sax_command);
                end
            end
        end
    end
end

fclose(fid_sax);
fclose(fid_sax_t2_map);
fclose(fid_sax_t2w);
