function apply_ai_cmr_landmark_detection_Perf(resDir, aiDir, pt_ids, files_record_picked, batch_size, sax_model, sax_model_pd, script_name, RO, E1, checkprocessed, suffix)
% apply_ai_cmr_landmark_detection_Perf(resDir, aiDir, pt_ids, files_record_picked, batch_size, sax_model, sax_model_pd, script_name, RO, E1, checkprocessed, suffix)

if(isunix())
    ext = '.sh';
else
    ext = '.bat';
end

fid_stress = fopen(fullfile(aiDir, [script_name '_cmr_landmark_stress' suffix ext]), 'a+');
fid_rest = fopen(fullfile(aiDir, [script_name '_cmr_landmark_rest' suffix ext]), 'a+');
fid_stress_pd = fopen(fullfile(aiDir, [script_name '_cmr_landmark_stress_pd' suffix ext]), 'a+');
fid_rest_pd = fopen(fullfile(aiDir, [script_name '_cmr_landmark_rest_pd' suffix ext]), 'a+');

set_env_apply_ai(fid_stress);
set_env_apply_ai(fid_rest);
set_env_apply_ai(fid_stress_pd);
set_env_apply_ai(fid_rest_pd);

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    if(isnumeric(pt_id))
        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' num2str(pt_id)]);
    else
        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    end
    
    case_stress = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'stress', 1);
    case_rest = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'rest', 1);
       
    if(numel(case_stress)>0)
        study_date = case_stress.study_dates(end,:);
        pt_id = case_stress.patientIDs{end,:};
    end
    if(numel(case_rest)>0)
        study_date = case_rest.study_dates(end,:);
        pt_id = case_rest.patientIDs{end,:};
    end
    disp([pt_id ' - ' study_date]);
    
    dst_dir = fullfile(aiDir, study_date, pt_id);
    if(~exist(dst_dir))
        dst_dir = fullfile(aiDir, pt_id);
    end
    
    % stress
    dst_dir_stress = fullfile(dst_dir, 'stress');
               
    % rest
    dst_dir_rest = fullfile(dst_dir, 'rest');
           
    contourDir = fullfile(dst_dir, 'res_ai');
       
    % lax_scripts
    if(isunix())
        command = ['python3 '];
    else
        command = ['python '];
    end
    
    if(exist(fullfile(dst_dir_stress, 'data_eigen.npy')))
        stress_command = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_stress, 'data_eigen.npy') ... 
            ' --output ' fullfile(contourDir, ['stress_AI_pts' suffix '.npy']) ' --prob ' fullfile(contourDir, ['stress_AI_probs' suffix '.npy']) ...
            ' --model ' sax_model ' --batch_size ' num2str(batch_size) ' --lax 0 --use_3D 0 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ];

        if(~checkprocessed || ~exist(fullfile(contourDir, ['stress_AI_pts' suffix '.npy'])))
            stress_command
            fprintf(fid_stress, '%s\n', stress_command);
        end
    end
    
    if(exist(fullfile(dst_dir_stress, 'pd.npy')))
        stress_command = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_stress, 'pd.npy') ... 
            ' --output ' fullfile(contourDir, ['stress_pd_AI_pts' suffix '.npy']) ' --prob ' fullfile(contourDir, ['stress_pd_AI_probs' suffix '.npy']) ...
            ' --model ' sax_model ' --batch_size ' num2str(batch_size) ' --lax 0 --use_3D 0 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ];

        if(~checkprocessed || ~exist(fullfile(contourDir, ['stress_pd_AI_pts' suffix '.npy'])))
            stress_command
            fprintf(fid_stress_pd, '%s\n', stress_command);
        end
    end
    
    if(exist(fullfile(dst_dir_rest, 'data_eigen.npy')))
        rest_command = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_rest, 'data_eigen.npy') ... 
            ' --output ' fullfile(contourDir, ['rest_AI_pts' suffix '.npy']) ' --prob ' fullfile(contourDir, ['rest_AI_probs' suffix '.npy']) ...
            ' --model ' sax_model ' --batch_size ' num2str(batch_size) ' --lax 0 --use_3D 0 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ];

        if(~checkprocessed || ~exist(fullfile(contourDir, ['rest_AI_pts' suffix '.npy'])))
            rest_command
            fprintf(fid_rest, '%s\n', rest_command);
        end
    end
    
    if(exist(fullfile(dst_dir_rest, 'pd.npy')))
        rest_command = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_rest, 'pd.npy') ... 
            ' --output ' fullfile(contourDir, ['rest_pd_AI_pts' suffix '.npy']) ' --prob ' fullfile(contourDir, ['rest_pd_AI_probs' suffix '.npy']) ...
            ' --model ' sax_model ' --batch_size ' num2str(batch_size) ' --lax 0 --use_3D 0 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ];

        if(~checkprocessed || ~exist(fullfile(contourDir, ['rest_pd_AI_pts' suffix '.npy'])))
            rest_command
            fprintf(fid_rest_pd, '%s\n', rest_command);
        end
    end
    
end

fclose(fid_stress);
fclose(fid_rest);
fclose(fid_stress_pd);
fclose(fid_rest_pd);
