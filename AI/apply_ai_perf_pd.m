function apply_ai_perf_pd(dataDir, resDir, aiDir, pt_ids, files_record_picked, model, script_name, RO, E1, checkprocessed, suffix)
% apply_ai_perf_pd(dataDir, resDir, aiDir, pt_ids, files_record_picked, model, script_name, RO, E1, checkprocessed, suffix)

if(isunix())
    ext = '.sh';
else
    ext = '.bat';
end

fid_stress_pd = fopen(fullfile(aiDir, [script_name '_cmr_perf_stress_pd' suffix ext]), 'a+');
fid_rest_pd = fopen(fullfile(aiDir, [script_name '_cmr_perf_rest_pd' suffix ext]), 'a+');

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
    
    if(exist(fullfile(dst_dir_stress, 'pd_for_seg.npy')))
        stress_command = [command 'perf_pd_seg.py --input ' fullfile(dst_dir_stress, 'pd_for_seg.npy') ... 
            ' --output ' fullfile(contourDir, ['stress_pd_seg_masks' suffix '.npy']) ' --prob ' fullfile(contourDir, ['stress_pd_seg_probs' suffix '.npy']) ...
            ' --contour ' fullfile(contourDir, ['stress_pd_seg_C' suffix '.npy']) ...
            ' --model ' model ' --cli_mode 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ' --thres 0.4 --contour_n_pts 64 --contour_smoothing 64' ];

        if(~checkprocessed || ~exist(fullfile(contourDir, ['stress_pd_seg_probs' suffix '.npy'])))
            stress_command
            fprintf(fid_stress_pd, '%s\n', stress_command);
        end
    end
    
    if(exist(fullfile(dst_dir_stress, 'pd_moco_for_seg.npy')))
        stress_command = [command 'perf_pd_seg.py --input ' fullfile(dst_dir_stress, 'pd_moco_for_seg.npy') ... 
            ' --output ' fullfile(contourDir, ['stress_pd_moco_seg_masks' suffix '.npy']) ' --prob ' fullfile(contourDir, ['stress_pd_moco_seg_probs' suffix '.npy']) ...
            ' --contour ' fullfile(contourDir, ['stress_pd_moco_seg_C' suffix '.npy']) ...
            ' --model ' model ' --cli_mode 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ' --thres 0.4 --contour_n_pts 64 --contour_smoothing 64' ];

        if(~checkprocessed || ~exist(fullfile(contourDir, ['stress_pd_moco_seg_probs' suffix '.npy'])))
            stress_command
            fprintf(fid_stress_pd, '%s\n', stress_command);
        end
    end
    
    if(exist(fullfile(dst_dir_rest, 'pd_for_seg.npy')))
        rest_command = [command 'perf_pd_seg.py --input ' fullfile(dst_dir_rest, 'pd_for_seg.npy') ... 
            ' --output ' fullfile(contourDir, ['rest_pd_seg_masks' suffix '.npy']) ' --prob ' fullfile(contourDir, ['rest_pd_seg_probs' suffix '.npy']) ...
            ' --contour ' fullfile(contourDir, ['rest_pd_seg_C' suffix '.npy']) ...
            ' --model ' model ' --cli_mode 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ' --thres 0.4 --contour_n_pts 64 --contour_smoothing 64' ];

        if(~checkprocessed || ~exist(fullfile(contourDir, ['rest_pd_seg_probs' suffix '.npy'])))
            rest_command
            fprintf(fid_rest_pd, '%s\n', rest_command);
        end
    end
    if(exist(fullfile(dst_dir_rest, 'pd_moco_for_seg.npy')))
        rest_command = [command 'perf_pd_seg.py --input ' fullfile(dst_dir_rest, 'pd_moco_for_seg.npy') ... 
            ' --output ' fullfile(contourDir, ['rest_pd_moco_seg_masks' suffix '.npy']) ' --prob ' fullfile(contourDir, ['rest_pd_moco_seg_probs' suffix '.npy']) ...
            ' --contour ' fullfile(contourDir, ['rest_pd_moco_seg_C' suffix '.npy']) ...
            ' --model ' model ' --cli_mode 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ' --thres 0.4 --contour_n_pts 64 --contour_smoothing 64' ];

        if(~checkprocessed || ~exist(fullfile(contourDir, ['rest_pd_moco_seg_probs' suffix '.npy'])))
            rest_command
            fprintf(fid_rest_pd, '%s\n', rest_command);
        end
    end
end

fclose(fid_stress_pd);
fclose(fid_rest_pd);

end
