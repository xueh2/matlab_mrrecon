function copy_results_ai_landmark_detection_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, dstDir, suffix, case_4chs, case_2chs, case_3chs, case_saxs)
% copy_results_ai_landmark_detection_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, dstDir, suffix)

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    
    if(numel(pt_ids)==size(case_4chs, 1))
        case_4ch = case_4chs{pt};
        case_2ch = case_2chs{pt};
        case_3ch = case_3chs{pt};
        case_sax = case_saxs{pt};
    else
        case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '4ch', 1);
        case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '2ch', 1);
        case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '3ch', 1);
        case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'sa', 1);
    end

    if(numel(case_4ch)==0 && numel(case_2ch)==0 && numel(case_3ch)==0)
        continue;
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

        src_dir = fullfile(aiDir, study_date, pt_id);
        if(~exist(src_dir))
            src_dir = fullfile(aiDir, pt_id);
            study_date = [];
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

        src_dir = fullfile(aiDir, pt_id);
        study_date = [];
    end
    contourDir = fullfile(src_dir, 'res_ai');
    
    dst_dir = fullfile(dstDir, study_date, pt_id);
    mkdir(dst_dir); 
    
    probs_str = ['_AI_probs' suffix '.npy'];
    pts_str = ['_AI_pts' suffix '.npy'];
    
    copy_one_view(case_4ch, resDir, src_dir, dst_dir, probs_str, pts_str, 'ch4', contourDir, suffix);
    copy_one_view(case_3ch, resDir, src_dir, dst_dir, probs_str, pts_str, 'ch3', contourDir, suffix);
    copy_one_view(case_2ch, resDir, src_dir, dst_dir, probs_str, pts_str, 'ch2', contourDir, suffix);
    
%     for d=1:size(case_4ch,1)
%         case_4ch_dir = fullfile(resDir, case_4ch.file_names{d});
%         [path, sname, ext] = fileparts(case_4ch_dir);
%         sname= sname(~isspace(sname));
%         src_dir_4ch = fullfile(src_dir, ['ch4_' sname]);
%         
%         if(exist(fullfile(src_dir_4ch, 'data.npy')) & exist(fullfile(contourDir, ['CH4' probs_str])))
%             dd = fullfile(dst_dir, ['ch4_' sname]);
%             mkdir(dd);
%             copyfile(fullfile(contourDir, ['CH4' probs_str]), dd);
%             copyfile(fullfile(contourDir, ['CH4' pts_str]), dd);
%             copyfile(fullfile(src_dir_4ch, 'data.npy'), dd);
%         end
%     end
    
%     
%     % 4ch
%     dst_dir_4ch = fullfile(dst_dir, 'ch4');
%                
%     % 2ch
%     dst_dir_2ch = fullfile(dst_dir, 'ch2');
%            
%     % 3ch
%     dst_dir_3ch = fullfile(dst_dir, 'ch3');
%     
%     %sax
%     dst_dir_sax = fullfile(dst_dir, 'sax');
%    
%     contourDir = fullfile(dst_dir, 'res_ai');
% 
%     dst_dir = fullfile(dstDir, study_date, pt_id);
%     mkdir(dst_dir); 
%     
%     probs_str = ['_AI_probs' suffix '.npy'];
%     pts_str = ['_AI_pts' suffix '.npy'];
%     
%     if(exist(fullfile(dst_dir_4ch, 'data.npy')) & exist(fullfile(contourDir, ['CH4' probs_str])))
%         dd = fullfile(dst_dir, 'ch4');
%         mkdir(dd);
%         copyfile(fullfile(contourDir, ['CH4' probs_str]), dd);
%         copyfile(fullfile(contourDir, ['CH4' pts_str]), dd);
%         copyfile(fullfile(dst_dir_4ch, 'data.npy'), dd);
%     end
%     
%     if(exist(fullfile(dst_dir_2ch, 'data.npy')) &exist(fullfile(contourDir, ['CH2' probs_str])))
%         dd = fullfile(dst_dir, 'ch2');
%         mkdir(dd);
%         copyfile(fullfile(contourDir, ['CH2' probs_str]), dd);
%         copyfile(fullfile(contourDir, ['CH2' pts_str]), dd);
%         copyfile(fullfile(dst_dir_2ch, 'data.npy'), dd);
%     end
%     
%     if(exist(fullfile(dst_dir_3ch, 'data.npy')) &exist(fullfile(contourDir, ['CH3' probs_str])))
%         dd = fullfile(dst_dir, 'ch3');
%         mkdir(dd);
%         copyfile(fullfile(contourDir, ['CH2' probs_str]), dd);
%         copyfile(fullfile(contourDir, ['CH3' pts_str]), dd);
%         copyfile(fullfile(dst_dir_3ch, 'data.npy'), dd);
%     end
%     
%     if(exist(fullfile(dst_dir_4ch, 'sax.npy')) &exist(fullfile(contourDir, 'SAX_AI_probs.npy')))
%         dd = fullfile(dst_dir, 'sax');
%         mkdir(dd);
%         copyfile(fullfile(contourDir, 'SAX_AI_probs.npy'), dd);
%         copyfile(fullfile(contourDir, 'SAX_AI_pts.npy'), dd);
%         copyfile(fullfile(dst_dir_sax, 'data.npy'), dd);
%     end
    
end
end

function copy_one_view(case_4ch, resDir, src_dir, dst_dir, probs_str, pts_str, view_str, contourDir, suffix)

    upper_view_str = upper(view_str);

    for d=1:size(case_4ch,1)
        case_4ch_dir = fullfile(resDir, case_4ch.file_names{d});
        [path, sname, ext] = fileparts(case_4ch_dir);
        sname= sname(~isspace(sname));
        src_dir_4ch = fullfile(src_dir, [view_str '_' sname]);
        if(~exist(fullfile(src_dir_4ch, 'data.npy')))
            src_dir_4ch = fullfile(src_dir, view_str);
        end
        
        pts_file = [];
        
        if(exist(fullfile(src_dir_4ch, 'data.npy')) & exist(fullfile(contourDir, [upper_view_str probs_str])))
            probs_file = fullfile(contourDir, [upper_view_str probs_str]);
            pts_file = fullfile(contourDir, [upper_view_str pts_str]);
        end
        
        if(exist(fullfile(src_dir_4ch, 'data.npy')) & exist(fullfile(contourDir, [upper_view_str '_AI_probs_' sname '_' suffix '.npy'])))            
            probs_file = fullfile(contourDir, [upper_view_str '_AI_probs_' sname '_' suffix '.npy']);
            pts_file = fullfile(contourDir, [upper_view_str '_AI_pts_' sname '_' suffix '.npy']);                        
        end
        
        if(~isempty(pts_file))
            pts = readNPY(pts_file);
            if(~isempty(find(pts_file(:)<0)))
                continue;
            end
            
            dd = fullfile(dst_dir, [view_str '_' sname]);
            mkdir(dd);
            copyfile(probs_file, fullfile(dd, [upper_view_str '_AI_probs.npy']));
            copyfile(pts_file, fullfile(dd, [upper_view_str '_AI_pts.npy']));
            copyfile(fullfile(src_dir_4ch, 'data.npy'), dd);
        end
    end
end