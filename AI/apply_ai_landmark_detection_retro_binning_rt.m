function apply_ai_landmark_detection_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, batch_size, lax_model, sax_model, three_D_model, script_name, RO, E1, RO_3D, E1_3D, checkprocessed, suffix)
% apply_ai_landmark_detection_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, batch_size, lax_model, sax_model, three_D_model, script_name, RO, E1, RO_3D, E1_3D, checkprocessed, suffix)

if(isunix())
    ext = '.sh';
else
    ext = '.bat';
end

fid_ch4 = fopen(fullfile(aiDir, [script_name '_cmr_landmark_ch4' ext]), 'a+');
fid_ch2 = fopen(fullfile(aiDir, [script_name '_cmr_landmark_ch2' ext]), 'a+');
fid_ch3 = fopen(fullfile(aiDir, [script_name '_cmr_landmark_ch3' ext]), 'a+');
fid_sax = fopen(fullfile(aiDir, [script_name '_cmr_landmark_sax' ext]), 'a+');

fid_ch4_3D = fopen(fullfile(aiDir, [script_name '_cmr_landmark_ch4_3D' ext]), 'a+');
fid_ch2_3D = fopen(fullfile(aiDir, [script_name '_cmr_landmark_ch2_3D' ext]), 'a+');
fid_ch3_3D = fopen(fullfile(aiDir, [script_name '_cmr_landmark_ch3_3D' ext]), 'a+');
fid_sax_3D = fopen(fullfile(aiDir, [script_name '_cmr_landmark_sax_3D' ext]), 'a+');

set_env_apply_ai(fid_ch4);
set_env_apply_ai(fid_ch4);
set_env_apply_ai(fid_ch2);
set_env_apply_ai(fid_ch3);
set_env_apply_ai(fid_sax);
set_env_apply_ai(fid_ch4_3D);
set_env_apply_ai(fid_ch2_3D);
set_env_apply_ai(fid_ch3_3D);
set_env_apply_ai(fid_sax_3D);

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};       
    
    if(isnumeric(pt_id))
        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' num2str(pt_id)]);
    else
        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    end
    
    if(numel(pt_ids)==size(case_4chs, 1))
        case_4ch = case_4chs{pt};
        case_2ch = case_2chs{pt};
        case_3ch = case_3chs{pt};
        case_sax = case_saxs{pt};
    else
        case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, {'4ch', 'ch4'});
        case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, {'2ch', 'ch2'});
        case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, {'3ch', 'ch3'});
    end
        
    try
        if(numel(case_4ch)>0)
            study_date = case_4ch.study_dates(end,:);
            pt_id = case_4ch.patientIDs{end,:};
        end
        if(numel(case_2ch)>0)
            study_date = case_2ch.study_dates(end,:);
            pt_id = case_2ch.patientIDs{end,:};
        end
        if(numel(case_3ch)>0)
            study_date = case_3ch.study_dates(end,:);
            pt_id = case_3ch.patientIDs{end,:};
        end

        disp([pt_id ' - ' study_date]);

        dst_dir = fullfile(aiDir, study_date, pt_id);
        if(~exist(dst_dir))
            dst_dir = fullfile(aiDir, pt_id);
        end    
    catch
        if(numel(case_4ch)>0)
            pt_id = case_4ch.patientIDs{end,:};
        end
        if(numel(case_2ch)>0)
            pt_id = case_2ch.patientIDs{end,:};
        end
        if(numel(case_3ch)>0)
            pt_id = case_3ch.patientIDs{end,:};
        end

        disp(pt_id);

        dst_dir = fullfile(aiDir, pt_id);
    end
           
    contourDir = fullfile(dst_dir, 'res_ai');
       
    % lax_scripts
    if(isunix())
        command = ['python3 '];
    else
        command = ['python '];
    end
    
    for d=1:size(case_4ch,1)
        case_4ch_dir = fullfile(resDir, case_4ch.file_names{d});
        [path, sname, ext] = fileparts(case_4ch_dir);
        sname= sname(~isspace(sname));
        dst_dir_4ch = fullfile(dst_dir, ['ch4_' sname]);
        if(~exist(fullfile(dst_dir_4ch, 'data.npy')))
            dst_dir_4ch = fullfile(dst_dir, ['ch4']);
        end
        if(exist(fullfile(dst_dir_4ch, 'data.npy')))
            ch4_command = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_4ch, 'data.npy') ... 
                ' --output ' fullfile(contourDir, ['CH4_AI_pts' '_' sname '_' suffix '.npy']) ' --prob ' fullfile(contourDir, ['CH4_AI_probs' '_' sname '_' suffix '.npy']) ...
                ' --model ' lax_model ' --batch_size ' num2str(batch_size) ' --lax 1 --use_3D 0 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ];

            ch4_command_3D = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_4ch, 'data.npy') ... 
                ' --output ' fullfile(contourDir, ['CH4_AI_pts_3D' '_' sname '_' suffix '.npy']) ...
                ' --prob ' fullfile(contourDir, ['CH4_AI_probs_3D' '_' sname '_' suffix '.npy']) ...
                ' --im_used ' fullfile(contourDir, ['CH4_AI_im_used_3D' '_' sname '_' suffix '.npy']) ...
                ' --model_3D ' three_D_model ' --lax 1 --use_3D 1 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO_3D) ' --E1 ' num2str(E1_3D) ];

            if(~checkprocessed || ~exist(fullfile(contourDir, ['CH4_AI_pts' '_' sname '_' suffix '.npy'])))
                ch4_command
                fprintf(fid_ch4, '%s\n', ch4_command);
            end

            if(~checkprocessed || ~exist(fullfile(contourDir, ['CH4_AI_pts_3D' '_' sname '_' suffix '.npy'])))
                ch4_command_3D
                fprintf(fid_ch4_3D, '%s\n', ch4_command_3D);
            end
        end
    end
    % -----------------------------
    
    for d=1:size(case_3ch,1)
        case_3ch_dir = fullfile(resDir, case_3ch.file_names{d});
        [path, sname, ext] = fileparts(case_3ch_dir);
        sname= sname(~isspace(sname));
        dst_dir_3ch = fullfile(dst_dir, ['ch3_' sname]);
        if(~exist(fullfile(dst_dir_3ch, 'data.npy')))
            dst_dir_3ch = fullfile(dst_dir, ['ch3']);
        end
        if(exist(fullfile(dst_dir_3ch, 'data.npy')))        
            ch3_command = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_3ch, 'data.npy') ... 
                ' --output ' fullfile(contourDir, ['CH3_AI_pts' '_' sname '_' suffix '.npy']) ' --prob ' fullfile(contourDir, ['CH3_AI_probs' '_' sname '_' suffix '.npy']) ...
                ' --model ' lax_model ' --batch_size ' num2str(batch_size) ' --lax 1 --use_3D 0 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ];

            ch3_command_3D = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_3ch, 'data.npy') ... 
                ' --output ' fullfile(contourDir, ['CH3_AI_pts_3D' '_' sname '_' suffix '.npy']) ...
                ' --prob ' fullfile(contourDir, ['CH3_AI_probs_3D' '_' sname '_' suffix '.npy']) ...
                ' --im_used ' fullfile(contourDir, ['CH3_AI_im_used_3D' '_' sname '_' suffix '.npy']) ...
                ' --model_3D ' three_D_model ' --lax 1 --use_3D 1 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO_3D) ' --E1 ' num2str(E1_3D) ];    

            if(~checkprocessed || ~exist(fullfile(contourDir, ['CH3_AI_pts' '_' sname '_' suffix '.npy'])))
                ch3_command
                fprintf(fid_ch3, '%s\n', ch3_command);
            end

            if(~checkprocessed || ~exist(fullfile(contourDir, ['CH3_AI_pts_3D' '_' sname '_' suffix '.npy'])))
                ch3_command_3D
                fprintf(fid_ch3_3D, '%s\n', ch3_command_3D);
            end
        end
    end
    % -----------------------------
    
    for d=1:size(case_2ch,1)
        case_2ch_dir = fullfile(resDir, case_2ch.file_names{d});
        [path, sname, ext] = fileparts(case_2ch_dir);
        sname= sname(~isspace(sname));
        dst_dir_2ch = fullfile(dst_dir, ['ch2_' sname]);
        if(~exist(fullfile(dst_dir_2ch, 'data.npy')))
            dst_dir_2ch = fullfile(dst_dir, ['ch2']);
        end
        if(exist(fullfile(dst_dir_2ch, 'data.npy')))   
            ch2_command = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_2ch, 'data.npy') ... 
                ' --output ' fullfile(contourDir, ['CH2_AI_pts' '_' sname '_' suffix '.npy']) ' --prob ' fullfile(contourDir, ['CH2_AI_probs' '_' sname '_' suffix '.npy']) ...
                ' --model ' lax_model ' --batch_size ' num2str(batch_size) ' --lax 1 --use_3D 0 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ];

            ch2_command_3D = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_2ch, 'data.npy') ... 
                ' --output ' fullfile(contourDir, ['CH2_AI_pts_3D' '_' sname '_' suffix '.npy']) ... 
                ' --prob ' fullfile(contourDir, ['CH2_AI_probs_3D' '_' sname '_' suffix '.npy']) ...
                ' --im_used ' fullfile(contourDir, ['CH2_AI_im_used_3D' '_' sname '_' suffix '.npy']) ...
                ' --model_3D ' three_D_model ' --lax 1 --use_3D 1 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO_3D) ' --E1 ' num2str(E1_3D) ];

            if(~checkprocessed || ~exist(fullfile(contourDir, ['CH2_AI_pts' '_' sname '_' suffix '.npy'])))
                ch2_command
                fprintf(fid_ch2, '%s\n', ch2_command);
            end

            if(~checkprocessed || ~exist(fullfile(contourDir, ['CH2_AI_pts_3D' '_' sname '_' suffix '.npy'])))
                ch2_command_3D
                fprintf(fid_ch2_3D, '%s\n', ch2_command_3D);
            end
        end
    end
    % -----------------------------
    
    % sax_scripts
%     if(exist(fullfile(dst_dir_sax, 'data.npy')))   
%         sax_command = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_sax, 'data.npy') ... 
%             ' --output ' fullfile(contourDir, ['SAX_AI_pts' '_' sname '_' suffix '.npy']) ' --prob ' fullfile(contourDir, ['SAX_AI_probs' '_' sname '_' suffix '.npy']) ...
%             ' --model ' sax_model ' --batch_size ' num2str(batch_size) ' --lax 0 --use_3D 0 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO) ' --E1 ' num2str(E1)];
% 
%         sax_command_3D = [command 'cmr_landmark_detection.py --input ' fullfile(dst_dir_sax, 'data.npy') ... 
%             ' --output ' fullfile(contourDir, ['SAX_AI_pts' '_' sname '_' suffix '.npy']) ' --prob ' fullfile(contourDir, ['SAX_AI_probs' '_' sname '_' suffix '.npy']) ...
%             ' --model_3D ' three_D_model ' --lax 0 --use_3D 1 --smooth_pts 0 --cli_mode 1 --fill_missed 1'];
%         
%         if(~checkprocessed || ~exist(fullfile(contourDir, ['SAX_AI_pts' '_' sname '_' suffix '.npy'])))
%             sax_command
%             fprintf(fid_sax, '%s\n', sax_command);
%         end
% 
%         if(~checkprocessed || ~exist(fullfile(contourDir, ['SAX_AI_pts_3D' '_' sname '_' suffix '.npy'])))
%             sax_command_3D
%             fprintf(fid_sax_3D, '%s\n', sax_command_3D);
%         end
%     end
    % -----------------------------
end

fclose(fid_ch4);
fclose(fid_ch3);
fclose(fid_ch2);
fclose(fid_sax);

fclose(fid_ch4_3D);
fclose(fid_ch3_3D);
fclose(fid_ch2_3D);
fclose(fid_sax_3D);
