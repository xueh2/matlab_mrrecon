function plot_results_ai_landmark_detection_T1_MOLLI(resDir, aiDir, pt_ids, files_record_picked)
% plot_results_ai_landmark_detection_T1_MOLLI(resDir, aiDir, pt_ids, files_record_picked)

    pic_dir = fullfile(aiDir, 'jpg_pics', 't1_molli_landmark_results');
    mkdir(pic_dir);

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
            
                if(exist(fullfile(contourDir, ['SAX_AI_pts_t1_map' '_plot.jpg'])))
                    continue;
                end
                
%                 if(exist(fullfile(t1_map_dir, 'data.npy')) & exist(view_file))
                  if(exist(fullfile(t1_map_dir, 'data.npy')))

%                     cmr_view = readNPY(view_file);                
%                     if(cmr_view==3) % SAX
                    
                        pts_file = fullfile(contourDir, ['SAX_AI_pts_t1_map.npy']);
                        probs_file = fullfile(contourDir, ['SAX_AI_probs_t1_map.npy']);
                        if(exist(pts_file))
                            data = readNPY(fullfile(t1_map_dir, 'data.npy'));
                            pts = readNPY(pts_file);
                            probs = readNPY(probs_file);                    
                            make_plots(data, pts, pic_dir, contourDir, case_used.file_names{k}, ['SAX_AI_pts_t1_map'], tt==2);
                        end
%                     end
                end            
            end
        end    
    end
end

function make_plots(data, pts, pic_dir, contourDir, data_file_name, prefix, is_post)
        
    h = figure;
    if(is_post) 
        WC = 550;
        WW = 650;        
    else
        WC = 1300;
        WW = 1300;
    end
    imagescn(data, [WC-WW/2 WC+WW/2], [], [12]);
    T1ColorMap;
    hold on
    pts = pts + 1;
    if(pts(1,1)>1)
        plot(pts(1,1), pts(1,2), 'r+', 'MarkerSize', 24, 'LineWidth', 4.0);
    end
    if(pts(2,1)>1)
        plot(pts(2,1), pts(2,2), 'b+', 'MarkerSize', 24, 'LineWidth', 4.0);
    end
    if(pts(3,1)>1)
        plot(pts(3,1), pts(3,2), 'w+', 'MarkerSize', 24, 'LineWidth', 4.0);
    end
    hold off
    
    saveas(h, fullfile(contourDir, [prefix '_plot.jpg']), 'jpg');
    saveas(h, fullfile(contourDir, [prefix '_plot.fig']), 'fig');
    
    copyfile(fullfile(contourDir, [prefix '_plot.jpg']), fullfile(pic_dir, [data_file_name '_' prefix '_plot.jpg']));
end

