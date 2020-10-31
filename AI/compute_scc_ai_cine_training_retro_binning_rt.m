function compute_scc_ai_cine_training_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_num, check_processed, plot_flag)
% compute_scc_ai_cine_training_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_num, check_processed, plot_flag)

mkdir(aiDir)

pic_dir = fullfile(aiDir, 'jpg_pics');
mkdir(pic_dir)
mkdir(fullfile(pic_dir, 'ch4_scc'))
mkdir(fullfile(pic_dir, 'ch2_scc'))
mkdir(fullfile(pic_dir, 'ch3_scc'))
mkdir(fullfile(pic_dir, 'sax_scc'))

visible_status = 'on';
if(~plot_flag)
    visible_status = 'off';
end

% pt_ids = unique(files_record_picked(:, 3));

ED = 1;
ES = 13;

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    
    case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '4ch');
    case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '2ch');
    case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '3ch');
    case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'sa');

    if(numel(case_4ch)==0 || numel(case_2ch)==0 || numel(case_sax)==0)
        continue;
    end
        
    case_prefix = [case_4ch.study_dates(end,:) '_' pt_id];
    dst_dir = fullfile(aiDir, case_4ch.study_dates(end,:), pt_id);
    
    ch4_processed = 0;
    ch2_processed = 0;
    ch3_processed = 0;
    sax_processed = 0;
    
    if(check_processed)
        dst_dir_sax = fullfile(dst_dir, 'sax');
        if(numel(case_3ch)==0)
            if(exist(fullfile(dst_dir_sax, 'record_header.mat')) ...
                & exist(fullfile(pic_dir, 'ch2_SCC', [case_prefix '_ch2_' 'ES.jpg'])) ...
                & exist(fullfile(pic_dir, 'ch4_SCC', [case_prefix '_ch4_' 'ES.jpg'])) ...
                & exist(fullfile(pic_dir, 'sax_SCC', [case_prefix '_sax_' 'ES_6.jpg'])) )
            
                disp(['already processed ' pt_id]);

                continue;
            end
        else
            if(exist(fullfile(dst_dir_sax, 'record_header.mat')) ...
                    & exist(fullfile(pic_dir, 'ch2_SCC', [case_prefix '_ch2_' 'ES.jpg'])) ...
                    & exist(fullfile(pic_dir, 'ch3_SCC', [case_prefix '_ch3_' 'ES.jpg'])) ...
                    & exist(fullfile(pic_dir, 'ch4_SCC', [case_prefix '_ch4_' 'ES.jpg'])) ...
                    & exist(fullfile(pic_dir, 'sax_SCC', [case_prefix '_sax_' 'ES_6.jpg'])) )

                disp(['already processed ' pt_id]);

                continue;
            end
        end
        
        if(exist(fullfile(pic_dir, 'ch4_SCC', [case_prefix '_ch4_' 'ES.npy'])))
            ch4_processed = 1;
        end
        
        if(exist(fullfile(pic_dir, 'ch2_SCC', [case_prefix '_ch2_' 'ES.npy'])))
            ch2_processed = 1;
        end
        
        if(exist(fullfile(pic_dir, 'ch3_SCC', [case_prefix '_ch3_' 'ES.npy'])))
            ch3_processed = 1;
        end
        
        if(exist(fullfile(pic_dir, 'sax_SCC', [case_prefix '_sax_' 'ED_6.npy'])))
            sax_processed = 1;
        end
    end
    
    try
        SLC = case_sax.headers{end}.encoding.encodingLimits.slice.maximum +1;
        sax_ind = size(case_sax, 1);
        if(SLC<6)
            SLC = case_sax.headers{1}.encoding.encodingLimits.slice.maximum +1;            
            sax_ind = 1;
        end

        % -----------------------------
        
        if(~ch4_processed)
            case_4ch_dir = fullfile(resDir, case_4ch.study_dates(end,:), case_4ch.file_names{end})
            [gt_4ch, gt_h_4ch, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_4ch_dir, series_num, 1);
            gt_4ch = squeeze(gt_4ch);
        end
        
        % -----------------------------
        
        if(~ch3_processed)
            if(numel(case_3ch)>0)
                case_3ch_dir = fullfile(resDir, case_3ch.study_dates(end,:), case_3ch.file_names{end})
                [gt_3ch, gt_h_3ch, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_3ch_dir, series_num, 1);
                gt_3ch = squeeze(gt_3ch);
            end
        end
        % -----------------------------
        
        if(~ch2_processed)
            case_2ch_dir = fullfile(resDir, case_2ch.study_dates(end,:), case_2ch.file_names{end})       
            [gt_2ch, gt_h_2ch, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_2ch_dir, series_num, 1);
            gt_2ch = squeeze(gt_2ch);
        end
        
        if(~sax_processed)
            case_sax_dir = fullfile(resDir, case_sax.study_dates(sax_ind,:), case_sax.file_names{sax_ind})       
            multi_slice_mode = 1;
            try   
                [gt_sax_all, gt_h_sax_all, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_sax_dir, series_num, 1);
                gt_sax_all = squeeze(gt_sax_all);
                multi_slice_mode = 0;
            catch
                rmdir(dst_dir, 's');
                continue;
            end

            if(size(gt_sax_all,4)==1 && SLC>1)
                multi_slice_mode = 1;
            end

%         if(numel(size(gt_sax_all))==4)
%             if(size(gt_sax_all,4)~= size(gt_2ch,3))
%                 rmdir(dst_dir, 's');
%                 continue;
%             end
%         end
        end
        
        if(plot_flag)
            h1 = figure; imagescn(permute(gt_4ch, [2 1 3]), [], [], [8], 3);
            h2 = figure; imagescn(permute(gt_2ch, [2 1 3]), [], [], [8], 3);
            h3 = figure; imagescn(permute(gt_3ch, [2 1 3]), [], [], [8], 3);            
            h4 = figure; imagescn(permute(gt_sax_all, [2 1 3 4]), [], [], [8], 4);
        end
    catch
        continue;
    end
    
    useMask = 1;
    boxFilterSize = 25;
    noisebackground = 30;
    thresRatioForNoise = 10;
    
    try
        mkdir(dst_dir)

        % 4ch
        if(~ch4_processed)
            dst_dir_4ch = fullfile(dst_dir, 'ch4');
            mkdir(dst_dir_4ch);

            if(numel(size(gt_4ch))==4)
                gt_4ch = squeeze(gt_4ch(:,:, round(size(gt_4ch,3)/2),:));
            end

            [a_eigen,V,D] = KL_Eigenimage(gt_4ch);
            [scc, mask] = Matlab_gt_surface_coil_correction(abs(double(a_eigen(:,:,end))), abs(double(a_eigen(:,:,end))), [], useMask, boxFilterSize, noisebackground, thresRatioForNoise);

            ind = find(scc==inf);
            scc(ind) = 1;
            b = gt_4ch ./ repmat(scc, [1 1 size(gt_4ch,3)]);
            b(ind) = 0;
            gt_4ch_norm = b/max(b(:)) * max(gt_4ch(:));

            if(plot_flag)
                figure; imagescn(permute(gt_4ch_norm, [2 1 3]), [], [], [8], 3);
            end 

            RO = size(gt_4ch,1)
            E1 = size(gt_4ch,2)
            ps = max(gt_h_4ch(1,1,1).FOV)/max(RO,E1)

            new_RO = round(RO*ps/dst_pixel_spacing(1));
            new_E1 = round(E1*ps/dst_pixel_spacing(2));

            data = Matlab_gt_resize_2D_image(double(gt_4ch), new_RO, new_E1, 5);
            data_dst = Matlab_gt_resize_2D_image(double(gt_4ch_norm), new_RO, new_E1, 5);
            scc_resampled = Matlab_gt_resize_2D_image(double(scc), new_RO, new_E1, 5);
            
            if(plot_flag)
                figure; imagescn(permute(data_dst, [2 1 3]), [], [], [8], 3);
            end

            gt_4ch = permute(gt_4ch, [2 1 3]);
            gt_4ch_norm = permute(gt_4ch_norm, [2 1 3]);

            data = permute(data, [2 1 3]);
            data_dst = permute(data_dst, [2 1 3]);
            size(data_dst)
            scc_resampled = permute(scc_resampled, [2 1]);

            RO = size(data_dst,1);
            E1 = size(data_dst,2);
            PHS = size(data_dst,3);

            data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), PHS);
            data_dst = zpad(data_dst,max(RO, NN_RO), max(E1, NN_E1), PHS);
            scc_resampled = zpad(scc_resampled,max(RO, NN_RO), max(E1, NN_E1));
            
            RO = size(data_dst,1);
            E1 = size(data_dst,2);
            PHS = size(data_dst,3);
            
            writeNPY(single(scc_resampled), fullfile(dst_dir_4ch, 'scc_resampled.npy'));
            writeNPY(single(scc_resampled), fullfile(pic_dir, 'ch4_SCC', [case_prefix '.npy']));
        end
        
        % 2ch
        if(~ch2_processed)
            dst_dir_2ch = fullfile(dst_dir, 'ch2');
            mkdir(dst_dir_2ch);

            if(numel(size(gt_2ch))==4)
                gt_2ch = squeeze(gt_2ch(:,:, round(size(gt_2ch,3)/2),:));
            end

            [a_eigen,V,D] = KL_Eigenimage(gt_2ch);
            [scc, mask] = Matlab_gt_surface_coil_correction(abs(double(a_eigen(:,:,end))), abs(double(a_eigen(:,:,end))), [], useMask, boxFilterSize, noisebackground, thresRatioForNoise);

            ind = find(scc==inf);
            scc(ind) = 1;
            b = gt_2ch ./ repmat(scc, [1 1 size(gt_2ch,3)]);
            b(ind) = 0;
            gt_2ch_norm = b/max(b(:)) * max(gt_2ch(:));

            if(plot_flag)
                figure; imagescn(permute(gt_2ch_norm, [2 1 3]), [], [], [8], 3);
            end

            RO = size(gt_2ch,1)
            E1 = size(gt_2ch,2)
            ps = max(gt_h_2ch(1,1,1).FOV)/max(RO,E1)

            new_RO = round(RO*ps/dst_pixel_spacing(1));
            new_E1 = round(E1*ps/dst_pixel_spacing(2));

            data = Matlab_gt_resize_2D_image(double(gt_2ch), new_RO, new_E1, 5);
            data_dst = Matlab_gt_resize_2D_image(double(gt_2ch_norm), new_RO, new_E1, 5);
            scc_resampled = Matlab_gt_resize_2D_image(double(scc), new_RO, new_E1, 5);

            if(plot_flag)
                figure; imagescn(permute(data_dst, [2 1 3]), [], [], [8], 3);
            end

            gt_2ch = permute(gt_2ch, [2 1 3]);
            gt_2ch_norm = permute(gt_2ch_norm, [2 1 3]);
            scc_resampled = permute(scc_resampled, [2 1]);

            data = permute(data, [2 1 3]);
            data_dst = permute(data_dst, [2 1 3]);
            size(data_dst)

            RO = size(data_dst,1);
            E1 = size(data_dst,2);
            PHS = size(data_dst,3);

            data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), PHS);
            data_dst = zpad(data_dst,max(RO, NN_RO), max(E1, NN_E1), PHS);
            scc_resampled = zpad(scc_resampled,max(RO, NN_RO), max(E1, NN_E1));

            RO = size(data_dst,1);
            E1 = size(data_dst,2);
            PHS = size(data_dst,3);

            if(PHS<ES)
                rmdir(dst_dir, 's');
                continue;
            end

            writeNPY(single(scc_resampled), fullfile(dst_dir_2ch, 'scc_resampled.npy'));
            writeNPY(single(scc_resampled), fullfile(pic_dir, 'ch2_SCC', [case_prefix '.npy']));           
        end
        
        % 3ch
        if(~ch3_processed && numel(case_3ch)>0)
            dst_dir_3ch = fullfile(dst_dir, 'ch3');
            mkdir(dst_dir_3ch);

            if(numel(size(gt_3ch))==4)
                gt_3ch = squeeze(gt_3ch(:,:, round(size(gt_3ch,3)/2),:));
            end

            useMask = 1;
            boxFilterSize = 25;
            noisebackground = 30;
            thresRatioForNoise = 10;

            [a_eigen,V,D] = KL_Eigenimage(gt_3ch);
            [scc, mask] = Matlab_gt_surface_coil_correction(abs(double(a_eigen(:,:,end))), abs(double(a_eigen(:,:,end))), [], useMask, boxFilterSize, noisebackground, thresRatioForNoise);

            ind = find(scc==inf);
            scc(ind) = 1;
            b = gt_3ch ./ repmat(scc, [1 1 size(gt_3ch,3)]);
            b(ind) = 0;
            gt_3ch_norm = b/max(b(:)) * max(gt_3ch(:));

            if(plot_flag)
                figure; imagescn(permute(gt_3ch_norm, [2 1 3]), [], [], [8], 3);
            end 

            RO = size(gt_3ch,1)
            E1 = size(gt_3ch,2)
            ps = max(gt_h_3ch(1,1,1).FOV)/max(RO,E1)

            new_RO = round(RO*ps/dst_pixel_spacing(1));
            new_E1 = round(E1*ps/dst_pixel_spacing(2));

            data = Matlab_gt_resize_2D_image(double(gt_3ch), new_RO, new_E1, 5);
            data_dst = Matlab_gt_resize_2D_image(double(gt_3ch_norm), new_RO, new_E1, 5);
            scc_resampled = Matlab_gt_resize_2D_image(double(scc), new_RO, new_E1, 5);

            if(plot_flag)
                figure; imagescn(permute(data_dst, [2 1 3]), [], [], [8], 3);
            end

            gt_3ch = permute(gt_3ch, [2 1 3]);
            gt_3ch_norm = permute(gt_3ch_norm, [2 1 3]);
            scc_resampled = permute(scc_resampled, [2 1]);

            data = permute(data, [2 1 3]);
            data_dst = permute(data_dst, [2 1 3]);
            size(data_dst)

            RO = size(data_dst,1);
            E1 = size(data_dst,2);
            PHS = size(data_dst,3);

            data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), PHS);
            data_dst = zpad(data_dst,max(RO, NN_RO), max(E1, NN_E1), PHS);
            scc_resampled = zpad(scc_resampled,max(RO, NN_RO), max(E1, NN_E1));

            RO = size(data_dst,1);
            E1 = size(data_dst,2);
            PHS = size(data_dst,3);

            if(PHS<ES)
                rmdir(dst_dir, 's');
                continue;
            end

            writeNPY(single(scc_resampled), fullfile(dst_dir_3ch, 'scc_resampled.npy'));
            writeNPY(single(scc_resampled), fullfile(pic_dir, 'ch3_SCC', [case_prefix '.npy']));
        end
    catch
        rmdir(dst_dir, 's');
        continue;
    end
        
    %sax
    scc_all = [];
    if(~sax_processed)
        dst_dir_sax = fullfile(dst_dir, 'sax');
        mkdir(dst_dir_sax);

        try
            gt_sax_norm = gt_sax_all;
            if(multi_slice_mode)
                gt_sax_norm = [];
                gt_sax_all = [];
                clear gt_h_sax_all
            end

            for slc=1:SLC

                disp(['process slice ' num2str(slc) ' ... ']);

                if(multi_slice_mode)
                    [gt_sax, gt_h_sax, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_sax_dir, (slc-1)*100, 1);
                    gt_sax = squeeze(gt_sax);

                    gt_h_sax = gt_h_sax(slc, :, :);

                    if(size(gt_sax,4)>1)
                        gt_sax = squeeze(gt_sax(:,:,slc,:));
                    end

                    gt_h_sax_all(slc) = gt_h_sax(1);
                else
                    gt_h_sax = gt_h_sax_all(:,slc);
                    gt_sax = squeeze(gt_sax_all(:,:,slc,:));
                end

                [a_eigen,V,D] = KL_Eigenimage(gt_sax);
                [scc, mask] = Matlab_gt_surface_coil_correction(abs(double(a_eigen(:,:,end))), abs(double(a_eigen(:,:,end))), [], useMask, boxFilterSize, noisebackground, thresRatioForNoise);
        
                ind = find(scc==inf);
                scc(ind) = 1;
                b = gt_sax ./ repmat(scc, [1 1 size(gt_sax,3)]);
                b(ind) = 0;
                gt_sax_norm(:,:,slc,:) = b/max(b(:)) * max(gt_sax(:));
                gt_sax_all(:,:,slc,:) = gt_sax;
                
                scc_all(:,:,slc) = scc;
            end

            if(plot_flag)
                figure; imagescn(permute(gt_sax_norm, [2 1 3 4]), [], [], [8], 4);
            end
        catch
            rmdir(dst_dir, 's');
            continue;
        end   

        % resample sax
        RO = size(gt_sax_all,1)
        E1 = size(gt_sax_all,2)
        ps = max(gt_h_sax_all(1,1,1).FOV)/max(RO,E1)

        new_RO = round(RO*ps/dst_pixel_spacing(1));
        new_E1 = round(E1*ps/dst_pixel_spacing(2));

        data = Matlab_gt_resize_2D_image(double(gt_sax_all), new_RO, new_E1, 5);
        data_dst = Matlab_gt_resize_2D_image(double(gt_sax_norm), new_RO, new_E1, 5);
        scc_all_resampled = Matlab_gt_resize_2D_image(double(scc_all), new_RO, new_E1, 5);

        if(plot_flag)
            figure; imagescn(permute(data_dst, [2 1 3 4]), [], [], [8], 4);
        end

        gt_sax_all = permute(gt_sax_all, [2 1 4 3]);
        gt_sax_norm = permute(gt_sax_norm, [2 1 4 3]);

        data = permute(data, [2 1 4 3]);
        size(data)
        data_dst = permute(data_dst, [2 1 4 3]);
        size(data_dst)
        scc_all_resampled = permute(scc_all_resampled, [2 1 3]);

        RO = size(data_dst,1);
        E1 = size(data_dst,2);
        PHS = size(data_dst,3);
        SLC = size(data_dst, 4);

        data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), PHS, SLC);
        data_dst = zpad(data_dst,max(RO, NN_RO), max(E1, NN_E1), PHS, SLC);
        scc_all_resampled = zpad(scc_all_resampled,max(RO, NN_RO), max(E1, NN_E1), SLC);

        RO = size(data_dst,1);
        E1 = size(data_dst,2);
        PHS = size(data_dst,3);
        SLC = size(data_dst, 4);

        if(PHS<ES)
            rmdir(dst_dir, 's');
            continue;
        end

        for slc=1:SLC
            writeNPY(single(scc_all_resampled(:,:,slc)), fullfile(pic_dir, 'sax_SCC', [case_prefix '_sax_' num2str(slc) '.npy']));
        end
    end
    
    closeall
    closeall    
end
