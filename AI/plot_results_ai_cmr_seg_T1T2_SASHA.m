function plot_results_ai_cmr_seg_T1T2_SASHA(resDir, aiDir, pt_ids, files_record_picked)
% plot_results_ai_cmr_seg_T1T2_SASHA(resDir, aiDir, pt_ids, files_record_picked)

    pic_dir = fullfile(aiDir, 'jpg_pics', 'cmr_seg_results');
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
            disp(['--> ' data_file_name]);
            contourDir = fullfile(dst_dir, [data_file_name '__res_ai']);

            t1w_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t1w_image']);
            t2w_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2w_image']);
            pd_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__pd_image']);
            t1_map_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t1_map']);
            t2_map_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__t2_map']);
            PSIR_map_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__PSIR_image']);
            DBPSIR_map_image_dir = fullfile(dst_dir, [case_used.file_names{k} '__DBPSIR_image']);
            
            if(~exist(t1w_image_dir))
                continue;
            end
            if(~exist(t2w_image_dir))
                continue;
            end
            if(~exist(pd_image_dir))
                continue;
            end
        
            endo_file = fullfile(contourDir, 'T1T2_seg_endo.npy');
            epi_file = fullfile(contourDir, 'T1T2_seg_epi.npy');
            if(~exist(endo_file))
                continue;
            end
            if(~exist(epi_file))
                continue;
            end
            
            t1_map = readNPY(fullfile(t1_map_image_dir, 'data.npy'));
            t2_map = readNPY(fullfile(t2_map_image_dir, 'data.npy'));
            t1w = readNPY(fullfile(t1w_image_dir, 'data.npy'));
            t2w = readNPY(fullfile(t2w_image_dir, 'data.npy'));
            pd = readNPY(fullfile(pd_image_dir, 'data.npy'));            
            PSIR = [];
            DBPSIR = [];
            if(exist(PSIR_map_image_dir))
                PSIR = readNPY(fullfile(PSIR_map_image_dir, 'data.npy'));
                record = load(fullfile(PSIR_map_image_dir, 'record.mat'));
                gt_h_PSIR = record.gt_h_view;
            end
            if(exist(DBPSIR_map_image_dir))
                DBPSIR = readNPY(fullfile(DBPSIR_map_image_dir, 'data.npy'));
                record = load(fullfile(DBPSIR_map_image_dir, 'record.mat'));
                gt_h_DBPSIR = record.gt_h_view;
            end
            
            endoC = readNPY(endo_file);
            epiC = readNPY(epi_file);
            
            mask = readNPY(fullfile(contourDir, 'T1T2_seg_mask.npy'));
            endo_mask = squeeze(mask(:,:,1,:));
            
            try
                pts = readNPY(fullfile(contourDir, 'SAX_AI_pts_t2w_image.npy'));
            catch
                pts = [];
            end
            
            SLC = size(t1_map, 3);
            
            t1v = [];
            for slc=1:SLC
                currMask = endo_mask(:,:,slc);
                currT1Map = t1_map(:,:,slc);
                ind = find(currMask>0);
                if(~isempty(ind))
                    v = currT1Map(ind);
                    t1v = [t1v; v];
                end
            end
            
            if(numel(t1v)>0)
                meanT1 = median(t1v(:));
            else
                meanT1 = 1300;
            end
            
            is_post = (meanT1<1000);
            
            make_plots(t1_map, endoC, epiC, pts, pic_dir, contourDir, case_used.file_names{k}, 't1_map', 1, 0, is_post, []);
            make_plots(t2_map, endoC, epiC, pts, pic_dir, contourDir, case_used.file_names{k}, 't2_map', 1, 1, is_post, []);
            make_plots(t1w, endoC, epiC, pts, pic_dir, contourDir, case_used.file_names{k}, 't1w', 0, 0, is_post, []);
            make_plots(t2w, endoC, epiC, pts, pic_dir, contourDir, case_used.file_names{k}, 't2w', 0, 0, is_post, []);
            make_plots(pd, endoC, epiC, pts, pic_dir, contourDir, case_used.file_names{k}, 'pd', 0, 0, is_post, []);
            
            if(~isempty(PSIR))
                make_plots(PSIR, endoC, epiC, pts, pic_dir, contourDir, case_used.file_names{k}, 'PSIR', 0, 0, 1, gt_h_PSIR);
            end
            if(~isempty(DBPSIR))
                make_plots(DBPSIR, endoC, epiC, pts, pic_dir, contourDir, case_used.file_names{k}, 'DBPSIR', 0, 0, 1, gt_h_DBPSIR);
            end
        end
    end
end

function make_plots(data, endoC, epiC, pts, pic_dir, contourDir, data_file_name, prefix, is_map, is_t2_map, is_post, gt_h)

    linewidth = 2.0;
    SLC = size(data, 3);

    h = figure;
    if(is_map) 
        if(is_post) 
            WC = 550;
            WW = 650;        
        else
            WC = 1300;
            WW = 1300;
        end
        
        if(is_t2_map)
            WC = 60;
            WW = 120;
        end
        imagescn(data, [WC-WW/2 WC+WW/2], [1 SLC], [18]);
        T1ColorMap;
    else
        v = median(data(:));
        WC = (0.5*v+7.5*v)/2;
        WW = (7.5*v-0.5*v);
        if(~isempty(gt_h))            
            if(gt_h(1,1).window_center>0)
                WC = gt_h(1,1).window_center;
                WW = gt_h(1,1).window_width;    
            end
        end
        imagescn(data, [WC-WW/2 WC+WW/2], [1 SLC], [18]);
    end
    
    h_axes=flipud(findobj(h,'type','axes'));
        
    for slc=1:SLC
        axes(h_axes(slc))
    
        if(~isempty(endoC))
            endo = endoC(:,:,slc);            
            if(endo(200,1)>0)
                endo = [endo; endo(1,:)];
                hold on  
                plot(endo(:,2)+1, endo(:,1)+1, 'y', 'LineWidth', linewidth);
                hold off
            end
        end
        
        if(~isempty(epiC))
            epi = epiC(:,:,slc);         
            if(epi(200,1)>0)
                epi = [epi; epi(1,:)];
                hold on    
                plot(epi(:,2)+1, epi(:,1)+1, 'g', 'LineWidth', linewidth);
                hold off
            end
        end
        
        hold on
        if(size(pts,3)==SLC)
            curr_pt = pts(:,:,slc);
            curr_pt = curr_pt + 1;
            if(curr_pt(1,1)>1)
                plot(curr_pt(1,1), curr_pt(1,2), 'r+', 'MarkerSize', 14, 'LineWidth', 2.0);
            end
            if(curr_pt(2,1)>1)
                plot(curr_pt(2,1), curr_pt(2,2), 'b+', 'MarkerSize', 14, 'LineWidth', 2.0);
            end
            if(curr_pt(3,1)>1)
                plot(curr_pt(3,1), curr_pt(3,2), 'k+', 'MarkerSize', 14, 'LineWidth', 2.0);
            end
        end        
        hold off
    end    
    
    saveas(h, fullfile(contourDir, [prefix '_plot.jpg']), 'jpg');
    saveas(h, fullfile(contourDir, [prefix '_plot.fig']), 'fig');
    
    copyfile(fullfile(contourDir, [prefix '_plot.jpg']), fullfile(pic_dir, [data_file_name '_' prefix '_plot.jpg']));
end

