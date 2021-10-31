function prepare_ai_cmr_training_T1T2_SASHA_ANTS(resDir, aiDir, pt_ids, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_t1num, series_t2num, series_num_t1map, series_num_t2map, series_t2num_moco, plot_flag)
% prepare_ai_cmr_training_T1T2_SASHA_ANTS(resDir, aiDir, pt_ids, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_t1num, series_t2num, series_num_t1map, series_num_t2map, series_t2num_moco, plot_flag)

mkdir(aiDir)

pic_dir = fullfile(aiDir, 'jpg_pics');
mkdir(pic_dir)
mkdir(fullfile(pic_dir, 't1_image'))
mkdir(fullfile(pic_dir, 't2_image'))
mkdir(fullfile(pic_dir, 't1_image_numpy'))
mkdir(fullfile(pic_dir, 't2_image_numpy'))
mkdir(fullfile(pic_dir, 't2_last_image'))
mkdir(fullfile(pic_dir, 't2_last_image_numpy'))
mkdir(fullfile(pic_dir, 't1_map'))
mkdir(fullfile(pic_dir, 't2_map'))
mkdir(fullfile(pic_dir, 't1_map_numpy'))
mkdir(fullfile(pic_dir, 't2_map_numpy'))

mkdir(fullfile(pic_dir, 'moco'))

mkdir(fullfile(pic_dir, 't1w_image'))
mkdir(fullfile(pic_dir, 't1w_image_numpy'))
mkdir(fullfile(pic_dir, 't2w_image'))
mkdir(fullfile(pic_dir, 't2w_image_numpy'))
mkdir(fullfile(pic_dir, 'pd_image'))
mkdir(fullfile(pic_dir, 'pd_image_numpy'))
mkdir(fullfile(pic_dir, 'PSIR_image'))
mkdir(fullfile(pic_dir, 'PSIR_image_numpy'))
mkdir(fullfile(pic_dir, 'DBPSIR_image'))
mkdir(fullfile(pic_dir, 'DBPSIR_image_numpy'))

visible_status = 'on';
if(~plot_flag)
    visible_status = 'off';
end

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
        
    case_prefix = [case_t1t2.study_dates(end,:) '_' pt_id];
    dst_dir = fullfile(aiDir, case_t1t2.study_dates(end,:), pt_id);
    
    pic_dir = fullfile(aiDir, 'jpg_pics');
    
    case_used = case_t1t2;
        
    for k=1:size(case_used, 1)
        dst_dir_case = fullfile(dst_dir, case_used.file_names{k});

        [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(case_used.file_names{k});
        
        if(str2num(measurementID)>300000)
            continue;
        end
        
        moco_fig = fullfile(pic_dir, 'moco', [case_used.file_names{k} '.fig']);
        t1_fig = fullfile(pic_dir, 'moco', [case_used.file_names{k} '_t1_map.fig']);
        t2_fig = fullfile(pic_dir, 'moco', [case_used.file_names{k} '_t2_map.fig']);
        if(exist(moco_fig) & exist(t1_fig) & exist(t2_fig))
            continue;
        end
        
        % process this case
        try
            case_dir = fullfile(resDir, case_used.study_dates(k,:), case_used.file_names{k})
            [gt_t1, gt_h_t1, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, series_t1num, 1);
            gt_t1 = squeeze(gt_t1);
            gt_h_t1 = squeeze(gt_h_t1);

            [gt_t1map, gt_h_t1map, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, series_num_t1map, 1);
            gt_t1map = squeeze(gt_t1map);
            gt_h_t1map = squeeze(gt_h_t1map);
           
            try
                [gt_t2, gt_h_t2, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, series_t2num, 1);
                gt_t2 = squeeze(gt_t2);
            catch
                [gt_t2, gt_h_t2, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, 101, 1);
                gt_t2 = squeeze(gt_t2);
            end
            gt_h_t2 = squeeze(gt_h_t2);
            
            [gt_t2map, gt_h_t2map, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, series_num_t2map, 1);
            gt_t2map = squeeze(gt_t2map);
            gt_h_t2map = squeeze(gt_h_t2map);
            
            [gt_t2_moco, gt_h_t2_moco, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, series_t2num_moco, 1);
            gt_t2_moco = squeeze(gt_t2_moco);
            gt_h_t2_moco = squeeze(gt_h_t2_moco);
                
            gt_PSIR = [];
            gt_h_PSIR = [];
            gt_DBPSIR = [];
            gt_DBPSIR = [];
                
%             try
%                 [gt_PSIR, gt_h_PSIR, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, series_PSIR, 1);
%                 gt_h_PSIR = squeeze(gt_h_PSIR);
%             catch
%                 gt_PSIR = [];
%                 gt_h_PSIR = [];
%             end
%             
%             try
%                 [gt_DBPSIR, gt_h_DBPSIR, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, series_DBPSIR, 1);
%                 gt_h_DBPSIR = squeeze(gt_h_DBPSIR);
%             catch
%                 gt_DBPSIR = [];
%                 gt_DBPSIR = [];
%             end
        catch
            continue;
        end

        if(max(gt_t1map(:))==0 | max(gt_t2map(:))==0)
            continue;
        end
        
%         if(numel(size(gt_t1))==5)
%             continue;
%         end
        
%         if(numel(size(gt_t1))==4 & numel(size(gt_t1map))==2)
%             gt_t1 = squeeze(gt_t1(:,:,1,:));
%             gt_t2 = squeeze(gt_t2(:,:,1,:));
%         end

        gt_t2_moco_ave = mean(gt_t2_moco, 5);
        t1w = gt_t2_moco(:,:,:,2);
        t2w = gt_t2_moco(:,:,:,end);
        pd = gt_t2_moco(:,:,:,1);

        SLC = size(gt_t2_moco, 3);
        
        %[t1w, t2w, pd] = compute_t1w_t2w_pd(gt_t2, gt_h_t2);

        gt_t1map = flipdim(gt_t1map, 2);
        gt_t2map = flipdim(gt_t2map, 2);
        t1w = flipdim(t1w, 2);
        t2w = flipdim(t2w, 2);
        pd = flipdim(pd, 2);
        gt_t2_moco_ave = flipdim(gt_t2_moco_ave, 2);
        
        for slc=1:SLC            
            if(size(gt_t1map,3)<slc)
                continue
            end

%                 prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't1_image'), 't1_image', squeeze(gt_t1(:,:,slc,:)), gt_h_t1, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, slc);
            prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't1_map'), 't1_map', squeeze(gt_t1map(:,:,slc)), gt_h_t1map, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, slc);

%                 prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't2_image'), 't2_image', squeeze(gt_t2(:,:,slc,:)), gt_h_t2, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, slc);
            prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't2_last_image'), 't2_last_image', squeeze(t2w(:,:,slc)), squeeze(gt_h_t2_moco(slc,1,1)), NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, slc);
            prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't2_map'), 't2_map', squeeze(gt_t2map(:,:,slc)), gt_h_t2map, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, slc);
        end

        prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't1_map'), 't1_map', gt_t1map, gt_h_t1map, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, -1);
        prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't2_map'), 't2_map', gt_t2map, gt_h_t2map, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, -1);
        
        prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't1w_image'), 't1w_image', t1w, gt_h_t2_moco(:,2,1), NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, -1);
        prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't2w_image'), 't2w_image', t2w, gt_h_t2_moco(:,end,1), NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, -1);
        prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 'pd_image'), 'pd_image', pd, gt_h_t2_moco(:,1,1), NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, -1);
        
        h = figure;
        imagescn(gt_t2_moco_ave, [0.1*mean(gt_t2_moco_ave(:)) 6*mean(gt_t2_moco_ave(:))], [1 SLC], [12], 4);        
        saveas(h, fullfile(pic_dir, 'moco', [case_used.file_names{k} '.fig']), 'fig');
        
        h = figure;
        imagescn(gt_t1map, [0 800], [1 SLC], [12]);        
        saveas(h, fullfile(pic_dir, 'moco', [case_used.file_names{k} '_t1_map.fig']), 'fig');
        
        h = figure;
        imagescn(gt_t2map, [0 120], [1 SLC], [12]); T1ColorMap;
        saveas(h, fullfile(pic_dir, 'moco', [case_used.file_names{k} '_t2_map.fig']), 'fig'); 

        closeall
        
%         if(~isempty(gt_PSIR))
%             prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 'PSIR_image'), 'PSIR_image', gt_PSIR, gt_h_PSIR, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, -1);
%         end
%         if(~isempty(gt_DBPSIR))
%             prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 'DBPSIR_image'), 'DBPSIR_image', gt_DBPSIR, gt_h_DBPSIR, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0, -1);
%         end
       
%         if(numel(size(gt_t1))==5)
%             gt_t1 = squeeze(gt_t1(:,:,1,:,1));
%             gt_t2 = squeeze(gt_t2(:,:,1,:,1));
%             gt_t1map = squeeze(gt_t1map(:,:,1,:));
%             gt_t2map = squeeze(gt_t2map(:,:,1,:));
%         end
%         
%         prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't1_image'), 't1_image', gt_t1, gt_h_t1, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0);
%         prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't1_map'), 't1_map', gt_t1map, gt_h_t1map, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0);
%         
%         prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't2_image'), 't2_image', gt_t2, gt_h_t2, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0);
%         prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't2_last_image'), 't2_last_image', gt_t2(:,:,end), gt_h_t2(end), NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0);
%         prepare_one_view(case_used.file_names{k}, dst_dir, fullfile(pic_dir, 't2_map'), 't2_map', gt_t2map, gt_h_t2map, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, 0);
    end
    
    closeall
    closeall    
end

end

function prepare_one_view(case_prefix, dst_dir, pic_dir, view_str, gt_view, gt_h_view, NN_RO, NN_E1, dst_pixel_spacing, plot_flag, visible_status, is_post, slc)

    if(nargin<13)
        slc = -1;
    end

    if(isempty(view_str))
        dst_dir_view = fullfile(dst_dir, case_prefix);
    else
        dst_dir_view = fullfile(dst_dir, [case_prefix '__' view_str]);
    end
    
    if(slc>0)
        dst_dir_view = [dst_dir_view '_slc' num2str(slc)];
    end
    
    if(~exist(dst_dir_view))
        mkdir(dst_dir_view);
    end
    
    RO = size(gt_view,1);
    E1 = size(gt_view,2);
    ps = max(gt_h_view(1,1).FOV)/max(RO,E1);

    WC = gt_h_view(1,1).window_center;
    WW = gt_h_view(1,1).window_width;

    if(WC<0)
        WC = mean(gt_view(:));
        WW =  mean(gt_view(:));
    end
    
    auto_windowing = 0;
    
    if(~isempty(view_str) & strcmp(view_str, 't1_map'))
        WC = 1300;
        WW = 1300;
    end
    
    if(~isempty(view_str) & (strcmp(view_str, 't2_last_image') | strcmp(view_str, 't1w_image') | strcmp(view_str, 't2w_image') | strcmp(view_str, 'pd_image')))
        WC = 2.5*median(gt_view(:));
        WW = 6*median(gt_view(:));
        auto_windowing = 1;
    end 
    
    if(~isempty(view_str) & strcmp(view_str, 't2_image'))
        WC = 2.5*median(gt_view(:));
        WW = 6*median(gt_view(:));
        auto_windowing = 1;
    end 
    
    if(~isempty(view_str) & strcmp(view_str, 't2_map'))
        WC = 60;
        WW = 120;
    end
    
    new_RO = round(RO*ps/dst_pixel_spacing(1));
    new_E1 = round(E1*ps/dst_pixel_spacing(2));

    if(mod(new_RO, 2)==1)
        new_RO = new_RO + 1;
    end
    
    if(mod(new_E1, 2)==1)
        new_E1 = new_E1 + 1;
    end
    
    gt_view(find(isnan(gt_view))) = 0;
    data = Matlab_gt_resize_2D_image(double(gt_view), new_RO, new_E1, 5);
    if(plot_flag)
        figure; imagescn(data, [WC-WW/2 WC+WW/2], [], []);
    end

    RO = size(data,1);
    E1 = size(data,2);
    SET = size(data,3);

    data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), SET);

%     data = permute(data, [2, 1, 3]);
   
    RO = size(data,1);
    E1 = size(data,2);
    SLC = size(data,3);

    data_resized_training = data;
    if(auto_windowing)
        v = 255 * data_resized_training/max(data_resized_training(:));
        for p=1:size(v,3)
            a2(:,:,p) = adapthisteq(uint8(v(:,:,p)));
        end
        h2 = figure('visible', visible_status); imagescn(a2);
    else
        h2 = figure('visible', visible_status); imagescn(data_resized_training, [WC-WW/2 WC+WW/2], [], []);
    end
    if(~isempty(view_str) & strcmp(view_str, 't1_map'))
        p = T1ColorMap;
        colormap(p);
    end
    if(~isempty(view_str) & strcmp(view_str, 't2_map'))
        p = T1ColorMap;
        colormap(p);
    end
    set(h2, 'visible', 'off')
    header = CreateGtImageHeader(data_resized_training);
    Matlab_gt_write_analyze(single(data_resized_training), header, fullfile(dst_dir_view, 'data_resized_training'));

    writeNPY(single(gt_view), fullfile(dst_dir_view, 'data_acq.npy'));
    writeNPY(single(data_resized_training), fullfile(dst_dir_view, 'data_resized_training.npy'));
    writeNPY(single(data), fullfile(dst_dir_view, 'data.npy'));
    saveas(h2, fullfile(dst_dir_view, 'data_resized_training'), 'jpg');

    pixel_spacing = ps;
    save(fullfile(dst_dir_view, 'record'), 'gt_h_view', 'dst_pixel_spacing', 'pixel_spacing', 'WC', 'WW');
    save(fullfile(dst_dir_view, 'record_header'), 'gt_h_view', 'dst_pixel_spacing', 'pixel_spacing', 'WC', 'WW');
    
    max_v = max(data_resized_training(:));
    
    for s=1:SET
        if(auto_windowing)
            pp = data_resized_training(:,:,s);
            a2 = adapthisteq(uint8(255 * pp/max_v));
            h2 = figure('visible', visible_status); imagescn(a2);
        else
            h2 = figure('visible', visible_status); imagescn(data_resized_training(:,:,s), [WC-WW/2 WC+WW/2], [], []);
        end
    
        if(~isempty(view_str) & strcmp(view_str, 't1_map'))
            p = T1ColorMap;
            colormap(p);
        end
        if(~isempty(view_str) & strcmp(view_str, 't2_map'))
            p = T1ColorMap;
            colormap(p);
        end
        if(slc>0)
            saveas(h2, fullfile(pic_dir, [case_prefix '_' num2str(s) '_slc' num2str(slc) '.jpg']), 'jpg');
            writeNPY(single(data_resized_training(:,:,s)), fullfile([pic_dir '_numpy'], [case_prefix '_' num2str(s) '_slc' num2str(slc) '.npy']));
        else            
            saveas(h2, fullfile(pic_dir, [case_prefix '_' num2str(s) '.jpg']), 'jpg');
            writeNPY(single(data_resized_training(:,:,s)), fullfile([pic_dir '_numpy'], [case_prefix '_' num2str(s) '.npy']));
        end
    end
    closeall
end

function [t1w, t2w, pd] = compute_t1w_t2w_pd(gt_t2, gt_h_t2)
    RO = size(gt_t2, 1);
    E1 = size(gt_t2, 2);
    
    if(numel(size(gt_t2))==4)
        SLC = size(gt_t2, 3);
        N = size(gt_t2, 4);
    else
        N = size(gt_t2, 3);
        gt_t2 = reshape(gt_t2, [RO, E1, 1, N]);
        SLC = 1;
        gt_h_t2 = reshape(gt_h_t2, [1, N]);
    end
    
    t1w = zeros(RO, E1, SLC);
    t2w = zeros(RO, E1, SLC);
    pd = zeros(RO, E1, SLC);
    
    for slc=1:SLC
        num_t1w = 0;
        num_t2w = 0;
        num_pd = 0;
        for n=1:N
            TE = gt_h_t2(slc, n).TE;
            TI = gt_h_t2(slc, n).TI;
            TS = gt_h_t2(slc, n).TS;
            if(TS>2000 & TE==0)
                num_pd = num_pd+1;
            end
            if(TS>2000 & TE==0)
                num_pd = num_pd+1;
                pd(:,:,slc) = pd(:,:,slc) + gt_t2(:,:,slc,n);
            end
            if(TS<=2000 & TE==0)
                num_t1w = num_t1w+1;
                t1w(:,:,slc) = t1w(:,:,slc) + gt_t2(:,:,slc,n);
            end
            if(TE>10)
                num_t2w = num_t2w+1;
                t2w(:,:,slc) = t2w(:,:,slc) + gt_t2(:,:,slc,n);
            end
        end
        
        if(num_t1w>0)
            t1w(:,:,slc) = t1w(:,:,slc) / num_t1w;
        end
        if(num_t2w>0)
            t2w(:,:,slc) = t2w(:,:,slc) / num_t2w;
        end
        if(num_pd>0)
            pd(:,:,slc) = pd(:,:,slc) / num_pd;
        end
    end
end