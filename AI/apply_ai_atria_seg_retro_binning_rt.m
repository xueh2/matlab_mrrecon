function apply_ai_atria_seg_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, ch4_model, script_name, RO, E1, checkprocessed, suffix)
% apply_ai_atria_seg_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, ch4_model, script_name, RO, E1, checkprocessed, suffix)

if(isunix())
    ext = '.sh';
else
    ext = '.bat';
end

fid_ch4 = fopen(fullfile(aiDir, [script_name '_cmr_atria_ch4' ext]), 'a+');
set_env_apply_ai(fid_ch4);

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
    
    script_name = 'cmr_atria_seg.py';
    
    for d=1:size(case_4ch,1)
        case_4ch_dir = fullfile(resDir, case_4ch.file_names{d});
        [path, sname, ext] = fileparts(case_4ch_dir);
        sname= sname(~isspace(sname));
        dst_dir_4ch = fullfile(dst_dir, ['ch4_' sname]);
        if(~exist(fullfile(dst_dir_4ch, 'data.npy')))
            dst_dir_4ch = fullfile(dst_dir, ['ch4']);
        end
        if(exist(fullfile(dst_dir_4ch, 'data.npy')))
            ch4_command = [command script_name ' --input ' fullfile(dst_dir_4ch, 'data.npy') ... 
                ' --prob ' fullfile(contourDir, ['CH4_AI_probs' '_' sname '_' suffix '.npy']) ...
                ' --C_LA ' fullfile(contourDir, ['CH4_AI_C_LA' '_' sname '_' suffix '.npy']) ...
                ' --C_RA ' fullfile(contourDir, ['CH4_AI_C_RA' '_' sname '_' suffix '.npy']) ...
                ' --fig ' fullfile(contourDir, ['CH4_AI_fig' '_' sname '_' suffix '.png']) ...
                ' --model ' ch4_model ' --cli_mode 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ];

            if(~checkprocessed || ~exist(fullfile(contourDir, ['CH4_AI_probs' '_' sname '_' suffix '.npy'])))
                ch4_command
                fprintf(fid_ch4, '%s\n', ch4_command);
            end
        end
    end
    % -----------------------------   
end

fclose(fid_ch4);
