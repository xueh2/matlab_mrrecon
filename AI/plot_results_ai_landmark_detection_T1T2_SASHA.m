function plot_results_ai_landmark_detection_T1T2_SASHA(resDir, aiDir, pt_ids, files_record_picked)
% plot_results_ai_landmark_detection_T1T2_SASHA(resDir, aiDir, pt_ids, files_record_picked)

    pic_dir = fullfile(aiDir, 'jpg_pics', 't1t2_landmark_results');
    mkdir(pic_dir);

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
            data_file_name = case_used.file_names{k};
            contourDir = fullfile(dst_dir, [data_file_name '__res_ai']);

            t2w_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2w_image']);        
            if(exist(fullfile(t2w_image_dir, 'data.npy')))
                pts_file = fullfile(contourDir, ['SAX_AI_pts_t2w_image.npy']);
                probs_file = fullfile(contourDir, ['SAX_AI_probs_t2w_image.npy']);
                if(exist(pts_file))
                    data = readNPY(fullfile(t2w_image_dir, 'data.npy'));
                    pts = readNPY(pts_file);
                    probs = readNPY(probs_file);                    
                    make_plots(data, pts, pic_dir, contourDir, data_file_name, ['SAX_AI_pts_t2w_image'], 0);
                end
            end
            
            t2_map_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2_map']);        
            if(exist(fullfile(t2_map_image_dir, 'data.npy')))
                pts_file = fullfile(contourDir, ['SAX_AI_pts_t2w_image.npy']);
                probs_file = fullfile(contourDir, ['SAX_AI_probs_t2w_image.npy']);
                if(exist(pts_file))
                    data = readNPY(fullfile(t2_map_image_dir, 'data.npy'));
                    pts = readNPY(pts_file);
                    probs = readNPY(probs_file);                    
                    make_plots(data, pts, pic_dir, contourDir, data_file_name, ['SAX_AI_pts_t2_map'], 1);
                end
            end
            
%             for slc=1:6
%                 t2_last_image_dir = fullfile(dst_dir, [data_file_name '__t2_last_image_slc' num2str(slc)]);
%                 t2_map_dir = fullfile(dst_dir, [data_file_name '__t2_map_slc' num2str(slc)]);
% 
%                 if(~exist(t2_last_image_dir))
%                     if(slc==1)
%                         t2_last_image_dir = fullfile(dst_dir, [data_file_name '__t2_last_image']);
%                         t2_map_dir = fullfile(dst_dir, [data_file_name '__t2_map']);
% 
%                         if(~exist(t2_last_image_dir))
%                             continue;
%                         end
%                     else
%                         continue;
%                     end
%                 end
%                 
%                 pts_file = fullfile(contourDir, ['SAX_AI_pts_t2_last_image_slc' num2str(slc) '.npy']);
%                 probs_file = fullfile(contourDir, ['SAX_AI_probs_t2_last_image_slc' num2str(slc) '.npy']);
%                 if(exist(pts_file))
%                     data = readNPY(fullfile(t2_last_image_dir, 'data.npy'));
%                     pts = readNPY(pts_file);
%                     probs = readNPY(probs_file);                    
%                     make_plots(data, pts, pic_dir, contourDir, data_file_name, ['SAX_AI_pts_t2_last_image_slc' num2str(slc)], 0);
%                 end
% 
%                 pts_file = fullfile(contourDir, ['SAX_AI_pts_t2_map_slc' num2str(slc) '.npy']);
%                 probs_file = fullfile(contourDir, ['SAX_AI_probs_t2_map_slc' num2str(slc) '.npy']);
%                 if(exist(pts_file))
%                     data = readNPY(fullfile(t2_map_dir, 'data.npy'));
%                     pts = readNPY(pts_file);
%                     probs = readNPY(probs_file);
%                     make_plots(data, pts, pic_dir, contourDir, data_file_name, ['SAX_AI_pts_t2_map_slc' num2str(slc)], 1);
%                 end
%                 
%                 closeall;
%             end
        end
    end
end

function make_plots(data, pts, pic_dir, contourDir, data_file_name, prefix, is_map)

    SLC = size(data, 3);

    h = figure;
    if(is_map)        
        imagescn(data, [0 120], [1 SLC], [18]);
        T1ColorMap;
    else
        v = median(data(:));
        imagescn(data, [0.5*v 7.5*v], [1 SLC], [18]);
    end
    
    h_axes=flipud(findobj(h,'type','axes'));
        
    for slc=1:SLC
        axes(h_axes(slc))
        hold on
        curr_pts = pts(:,:,slc) + 1;
        if(curr_pts(1,1)>1)
            plot(curr_pts(1,1), curr_pts(1,2), 'r+', 'MarkerSize', 16, 'LineWidth', 2.0);
        end
        if(curr_pts(2,1)>1)
            plot(curr_pts(2,1), curr_pts(2,2), 'b+', 'MarkerSize', 16, 'LineWidth', 2.0);
        end
        if(curr_pts(3,1)>1)
            plot(curr_pts(3,1), curr_pts(3,2), 'k+', 'MarkerSize', 16, 'LineWidth', 2.0);
        end
        hold off
    end
    
    saveas(h, fullfile(contourDir, [prefix '_plot.jpg']), 'jpg');
    saveas(h, fullfile(contourDir, [prefix '_plot.fig']), 'fig');
    
    copyfile(fullfile(contourDir, [prefix '_plot.jpg']), fullfile(pic_dir, [data_file_name '_' prefix '_plot.jpg']));
end

