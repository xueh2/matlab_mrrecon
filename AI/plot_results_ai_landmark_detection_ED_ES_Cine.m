function plot_results_ai_landmark_detection_ED_ES_Cine(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, check_processed, suffix)
% plot_results_ai_landmark_detection_ED_ES_Cine(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, check_processed, suffix)

    pic_dir = fullfile(aiDir, 'jpg_pics', 'Cine_landmark_results');
    mkdir(pic_dir);

    ED = 1;
    ES = 13;

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

        if(numel(case_4ch)==0 & numel(case_3ch)==0 & numel(case_2ch)==0)
            continue;
        end

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

        % 4ch
        dst_dir_4ch = fullfile(dst_dir, 'ch4');

        % 2ch
        dst_dir_2ch = fullfile(dst_dir, 'ch2');

        % 3ch
        dst_dir_3ch = fullfile(dst_dir, 'ch3');

        %sax
        dst_dir_sax = fullfile(dst_dir, 'sax');

        contourDir = fullfile(dst_dir, 'res_ai');
       
        for v=1:4
            if(v==1)
                if(numel(case_4ch)==0)
                    continue;
                end
                dst_dir_view = dst_dir_4ch;
                prefix = 'CH4';
                case_used = case_4ch(end,:);
            end
            if(v==2)                
                if(numel(case_2ch)==0)
                    continue;
                end
                dst_dir_view = dst_dir_2ch;
                prefix = 'CH2';
                case_used = case_2ch(end,:);
            end
            if(v==3)                
                if(numel(case_3ch)==0)
                    continue;
                end
                dst_dir_view = dst_dir_3ch;
                prefix = 'CH3';
                case_used = case_3ch(end,:);
            end
            if(v==4)                
                if(numel(case_sax)==0)
                    continue;
                end
                dst_dir_view = dst_dir_sax;
                prefix = 'SAX';
                case_used = case_sax(end,:);
            end
            
            if(~exist(fullfile(dst_dir_view, 'data.npy')))
                continue;
            end
            
            if(check_processed & exist(fullfile(contourDir, [prefix '_plot_ES' suffix '.jpg'])))
                continue;
            end
            
            record = load(fullfile(dst_dir_view, 'record'));
            pts_file = fullfile(contourDir, [prefix '_AI_pts' suffix '.npy']);
            
            if(~exist(pts_file))
                continue;
            end
            
            probs_file = fullfile(contourDir, [prefix  '_AI_probs' suffix '.npy']);
            data = readNPY(fullfile(dst_dir_view, 'data.npy'));
            pts = readNPY(pts_file);
            probs = readNPY(probs_file);                    
            make_plots(squeeze(data(:,:,ED,:)), squeeze(pts(:,:,ED,:)), pic_dir, contourDir, case_used.file_names{1}, prefix, ['ED' suffix]);
            make_plots(squeeze(data(:,:,ES,:)), squeeze(pts(:,:,ES,:)), pic_dir, contourDir, case_used.file_names{1}, prefix, ['ES' suffix]);
        end
        
        closeall
    end
end

function make_plots(data, pts, pic_dir, contourDir, data_file_name, prefix, phs_str)

    SLC = size(data, 3);

    if(SLC==1)
        h = figure;
        imagescn(data/max(data(:)), [0 0.45], [], [8]);
        hold on
        pts = pts + 1;
        if(pts(1,1)>1)
            plot(pts(1,1), pts(1,2), 'r.', 'MarkerSize', 24, 'LineWidth', 2.0);
        end
        if(pts(2,1)>1)
            plot(pts(2,1), pts(2,2), 'b.', 'MarkerSize', 24, 'LineWidth', 2.0);
        end
        if(pts(3,1)>1)
            plot(pts(3,1), pts(3,2), 'w.', 'MarkerSize', 24, 'LineWidth', 2.0);
        end
        hold off
    else
        if(mod(SLC,2)==0)
            h = figure;imagescn(data/max(data(:)), [0 0.45], [2 SLC/2], [24]);
        else
            h = figure;imagescn(data/max(data(:)), [0 0.45], [2 floor(SLC/2)+1], [24]);
        end
        
        h_axes = findobj(h,'type','axes');
        
        for slc=1:SLC
            axes(h_axes(slc));
            curr_pts = pts(:,:,SLC-slc+1);
            hold on
            curr_pts = curr_pts + 1;
            if(curr_pts(1,1)>1)
                plot(curr_pts(1,1), curr_pts(1,2), 'r+', 'MarkerSize', 12, 'LineWidth', 2.0);
            end
            if(curr_pts(2,1)>1)
                plot(curr_pts(2,1), curr_pts(2,2), 'b+', 'MarkerSize', 12, 'LineWidth', 2.0);
            end
            if(curr_pts(3,1)>1)
                plot(curr_pts(3,1), curr_pts(3,2), 'w+', 'MarkerSize', 12, 'LineWidth', 2.0);
            end
            hold off
        end
    end
    
    saveas(h, fullfile(contourDir, [prefix '_plot_' phs_str '.jpg']), 'jpg');
    saveas(h, fullfile(contourDir, [prefix '_plot_' phs_str '.fig']), 'fig');
    
    copyfile(fullfile(contourDir, [prefix '_plot_' phs_str '.jpg']), fullfile(pic_dir, [data_file_name '_' prefix '_plot_' phs_str '.jpg']));
end

