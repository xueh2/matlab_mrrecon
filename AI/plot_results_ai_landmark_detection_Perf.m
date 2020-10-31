function plot_results_ai_landmark_detection_Perf(resDir, aiDir, pt_ids, files_record_picked, suffix)
% plot_results_ai_landmark_detection_Perf(resDir, aiDir, pt_ids, files_record_picked, suffix)

    pic_dir = fullfile(aiDir, 'jpg_pics', 'perf_eigen_landmark_results');
    pic_pd_dir = fullfile(aiDir, 'jpg_pics', 'perf_pd_landmark_results');
    mkdir(pic_dir);
    mkdir(pic_pd_dir);

    for pt=1:numel(pt_ids)
    
        closeall

        pt_id = pt_ids{pt};

        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);

        case_stress = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'stress');
        case_rest = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'rest');

        if(case_stress.headers{1}.subjectInformation.patientGender == 'O')
            continue;
        end

        dst_dir = fullfile(aiDir, case_stress.study_dates(end,:), pt_id);
        contourDir = fullfile(dst_dir, 'res_ai');

        dst_dir_stress = fullfile(dst_dir, 'stress');
        dst_dir_rest = fullfile(dst_dir, 'rest');
    
        pts_file = fullfile(contourDir, ['stress_AI_pts' suffix '.npy']);
        probs_file = fullfile(contourDir, ['stress_AI_probs' suffix '.npy']);
        if(exist(pts_file))
            data = readNPY(fullfile(dst_dir_stress, 'data_eigen.npy'));
            pts = readNPY(pts_file);
            probs = readNPY(probs_file);
            SLC = size(data, 3);
            for slc=1:SLC
                make_plots(data(:,:,slc), pts(:,:,slc), pic_dir, contourDir, case_stress.file_names{end,:}, ['sterss_AI_pts_eigen_slc' num2str(slc)], 0, 0);
            end
        end
        
        pts_file = fullfile(contourDir, ['rest_AI_pts' suffix '.npy']);
        probs_file = fullfile(contourDir, ['rest_AI_probs' suffix '.npy']);
        if(exist(pts_file))
            data = readNPY(fullfile(dst_dir_rest, 'data_eigen.npy'));
            pts = readNPY(pts_file);
            probs = readNPY(probs_file);
            SLC = size(data, 3);
            for slc=1:SLC
                make_plots(data(:,:,slc), pts(:,:,slc), pic_dir, contourDir, case_rest.file_names{end,:}, ['rest_AI_pts_eigen_slc' num2str(slc)], 0, 0);
            end
        end
        
        % ---------------------------------------
        
        pts_file = fullfile(contourDir, ['stress_pd_AI_pts' suffix '.npy']);
        probs_file = fullfile(contourDir, ['stress_pd_AI_probs' suffix '.npy']);
        if(exist(pts_file))
            data = readNPY(fullfile(dst_dir_stress, 'pd.npy'));
            pts = readNPY(pts_file);
            probs = readNPY(probs_file);
            SLC = size(data, 3);
            for slc=1:SLC
                make_plots(data(:,:,slc), pts(:,:,slc), pic_pd_dir, contourDir, case_stress.file_names{end,:}, ['sterss_AI_pts_PD_slc' num2str(slc)], 0, 1);
            end
        end
        
        pts_file = fullfile(contourDir, ['rest_pd_AI_pts' suffix '.npy']);
        probs_file = fullfile(contourDir, ['rest_pd_AI_probs' suffix '.npy']);
        if(exist(pts_file))
            data = readNPY(fullfile(dst_dir_rest, 'pd.npy'));
            pts = readNPY(pts_file);
            probs = readNPY(probs_file);
            SLC = size(data, 3);
            for slc=1:SLC
                make_plots(data(:,:,slc), pts(:,:,slc), pic_pd_dir, contourDir, case_rest.file_names{end,:}, ['rest_AI_pts_PD_slc' num2str(slc)], 0, 1);
            end
        end
    end
end

function make_plots(data, pts, pic_dir, contourDir, data_file_name, prefix, is_map, only_C_LV)

    data = squeeze(data);
    pts = squeeze(pts);

    h = figure;
    if(is_map)        
        imagescn(data, [0 400], [], [8]);
        T1ColorMap;
    else
        v = median(data(:));
        imagescn(data, [0.25*v 10*v], [], [8]);
    end
    hold on
    pts = pts + 1;
    if(~only_C_LV & pts(1,1)>1)
        plot(pts(1,1), pts(1,2), 'r+', 'MarkerSize', 24, 'LineWidth', 4.0);
    end
    if(~only_C_LV & pts(2,1)>1)
        plot(pts(2,1), pts(2,2), 'b+', 'MarkerSize', 24, 'LineWidth', 4.0);
    end
    if(pts(3,1)>1)
        plot(pts(3,1), pts(3,2), 'y+', 'MarkerSize', 24, 'LineWidth', 4.0);
    end
    hold off
    
    saveas(h, fullfile(contourDir, [prefix '_plot.jpg']), 'jpg');
    saveas(h, fullfile(contourDir, [prefix '_plot.fig']), 'fig');
    
    copyfile(fullfile(contourDir, [prefix '_plot.jpg']), fullfile(pic_dir, [data_file_name '_' prefix '_plot.jpg']));
end

