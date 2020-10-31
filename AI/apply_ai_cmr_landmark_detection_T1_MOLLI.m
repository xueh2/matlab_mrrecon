function apply_ai_cmr_landmark_detection_T1_MOLLI(resDir, aiDir, pt_ids, files_record_picked, sax_model, script_name, checkprocessed)
% apply_ai_cmr_landmark_detection_T1_MOLLI(resDir, aiDir, pt_ids, files_record_picked, sax_model, script_name, checkprocessed)

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
            
fid_sax = fopen(fullfile(aiDir, [script_name '_cmr_landmark_sax' ext]), 'a+');

set_env_apply_ai(fid_sax);

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
            contourDir = fullfile(dst_dir, [case_used.file_names{k} '__res_ai']);
            view_file = fullfile(contourDir, ['CMR_view_res_t1_map.npy']);
            
            %if(exist(fullfile(t1_map_dir, 'data.npy')) & exist(view_file))
            if(exist(fullfile(t1_map_dir, 'data.npy')))
                                
%                 cmr_view = readNPY(view_file);
                
%                 if(cmr_view==3) % SAX
                    pts_file = fullfile(contourDir, ['SAX_AI_pts_t1_map.npy']);
                    probs_file = fullfile(contourDir, ['SAX_AI_probs_t1_map.npy']);
                    sax_command = [command_pt ' cmr_landmark_detection.py --input ' fullfile(t1_map_dir, 'data.npy') ... 
                        ' --output ' pts_file ' --prob ' probs_file ...
                        ' --model ' sax_model ' --batch_size 16' ' --lax 0 --use_3D 0 --cli_mode 1 --RO 352 --E1 352 --smooth_pts 0'];

                    if(~checkprocessed || ~exist(pts_file))
                        sax_command
                        fprintf(fid_sax, '%s\n', sax_command);
                    end
%                 end
            end            
        end
    end
end

fclose(fid_sax);
