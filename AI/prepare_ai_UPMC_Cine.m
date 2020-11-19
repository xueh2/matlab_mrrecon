function prepare_ai_UPMC_Cine(dataDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, dst_pixel_spacing, NN_RO, NN_E1, use_moco, regularization_hilbert_strength, check_processed, plot_flag)
% prepare_ai_UPMC_Cine(dataDir, aiDir, pt_ids, files_record_picked, dst_pixel_spacing, NN_RO, NN_E1, series_num, use_moco, regularization_hilbert_strength, check_processed, plot_flag)

    mkdir(aiDir)

    pic_dir = fullfile(aiDir, 'jpg_pics');
    mkdir(pic_dir)
    mkdir(fullfile(pic_dir, 'ch4'))
    mkdir(fullfile(pic_dir, 'ch2'))
    mkdir(fullfile(pic_dir, 'ch3'))
    mkdir(fullfile(pic_dir, 'sax'))
    mkdir(fullfile(pic_dir, 'ch4_original'))
    mkdir(fullfile(pic_dir, 'ch2_original'))
    mkdir(fullfile(pic_dir, 'ch3_original'))
    mkdir(fullfile(pic_dir, 'sax_original'))
    mkdir(fullfile(pic_dir, 'ch4_numpy'))
    mkdir(fullfile(pic_dir, 'ch2_numpy'))
    mkdir(fullfile(pic_dir, 'ch3_numpy'))
    mkdir(fullfile(pic_dir, 'sax_numpy'))
    mkdir(fullfile(pic_dir, 'ch4_original_numpy'))
    mkdir(fullfile(pic_dir, 'ch2_original_numpy'))
    mkdir(fullfile(pic_dir, 'ch3_original_numpy'))
    mkdir(fullfile(pic_dir, 'sax_original_numpy'))

    visible_status = 'on';
    if(~plot_flag)
        visible_status = 'off';
    end

    % pt_ids = unique(files_record_picked(:, 3));

    phases = [1:10:30, 13, 30];

    for pt=1:numel(pt_ids)

        closeall

        pt_id = pt_ids{pt};

        disp(['-------> ' num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);

        if(numel(pt_ids)==size(case_4chs, 1))
            case_4ch = case_4chs{pt};
            case_2ch = case_2chs{pt};
            case_3ch = case_3chs{pt};
            case_sax = case_saxs{pt};
        else
            case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'CH4', 0);
            case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'CH2', 0);
            case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'CH3', 0);
            case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'SAX', 0);
        end
        
        if(numel(case_4ch)==0 & numel(case_2ch)==0 & numel(case_3ch)==0)
            continue;
        end

%         if(numel(case_4ch)>0)
%             study_date = case_4ch.study_dates(end,:);
%         end
%         if(numel(case_2ch)>0)
%             study_date = case_2ch.study_dates(end,:);
%         end
%         if(numel(case_3ch)>0)
%             study_date = case_3ch.study_dates(end,:);
%         end

%         pt_id = str2num(pt_id);
%         pt_id = sprintf('%04d', pt_id);

        case_prefix = [pt_id]
        dst_dir = fullfile(aiDir, pt_id);

        ch4_processed = 0;
        ch2_processed = 0;
        ch3_processed = 0;
        sax_processed = 0;

        try
            % -----------------------------

            dst_dir = fullfile(aiDir, pt_id);
            
            if(~ch4_processed && numel(case_4ch)>0)
                process_one_view_all(dataDir, dst_dir, case_4ch, check_processed, pt_id, pic_dir, 'ch4', dst_pixel_spacing, NN_RO, NN_E1, use_moco, regularization_hilbert_strength, phases, plot_flag, visible_status);
            end

            % -----------------------------

            if(~ch3_processed && numel(case_3ch)>0)
                process_one_view_all(dataDir, dst_dir, case_3ch, check_processed, pt_id, pic_dir, 'ch3', dst_pixel_spacing, NN_RO, NN_E1, use_moco, regularization_hilbert_strength, phases, plot_flag, visible_status);
            end
            % -----------------------------

            if(~ch2_processed && numel(case_2ch)>0)
                process_one_view_all(dataDir, dst_dir, case_2ch, check_processed, pt_id, pic_dir, 'ch2', dst_pixel_spacing, NN_RO, NN_E1, use_moco, regularization_hilbert_strength, phases, plot_flag, visible_status);
            end
        catch
            continue;
        end

        closeall
        closeall    
    end

end

function process_one_view_all(dataDir, dst_dir, case_4ch, check_processed, pt_id, pic_dir, view_str, dst_pixel_spacing, NN_RO, NN_E1, use_moco, regularization_hilbert_strength, phases, plot_flag, visible_status)
    for d=1:size(case_4ch,1)
        case_4ch_dir = fullfile(dataDir, case_4ch.file_names{d});
        if(exist(case_4ch_dir))
            disp([ '---> ' case_4ch.file_names{d}]);
            [path, sname, ext] = fileparts(case_4ch_dir); 
            sname= sname(~isspace(sname));
            case_prefix = [view_str '_' sname];
            if(check_processed)
                if(exist(fullfile(dst_dir, case_prefix, 'data.npy')))
                    continue;
                end
            end                        

            gt_4ch = readNPY(fullfile(case_4ch_dir, 'data.npy'));
            gt_4ch_ps = readNPY(fullfile(case_4ch_dir, 'PixelSpacing.npy'));
            gt_4ch = squeeze(gt_4ch);
            if(size(gt_4ch, 4)>1)
                SLC = size(gt_4ch, 4);
                figure; imagescn(gt_4ch, [], [1 SLC], [], 3);
                reply = input('Which slice to pick :','s');
                if isempty(reply)
                    slc = round(SLC/2);
                else
                    slc = str2num(reply);
                end
                gt_4ch = gt_4ch(:,:,:, slc);
                closeall
            end

            process_one_view([pt_id '_' case_prefix], dst_dir, fullfile(pic_dir, view_str), fullfile(pic_dir, [view_str '_numpy']), [view_str '_' sname], gt_4ch, gt_4ch_ps, dst_pixel_spacing, NN_RO, NN_E1, use_moco, regularization_hilbert_strength, phases, plot_flag, visible_status);
%             save(fullfile(dst_dir, ['case_' view_str '_dir_' sname '.mat']), ['case_' view_str '_' sname '_dir']);

            if(plot_flag)
                h1 = figure; imagescn(gt_4ch, [], [], [8], 3);
            end
        end
    end
end

function process_one_view(case_prefix, dst_dir, pic_dir, pic_dir_numpy, view_str, gt_4ch, gt_4ch_ps, dst_pixel_spacing, NN_RO, NN_E1, use_moco, regularization_hilbert_strength, phases, plot_flag, visible_status)
    try        
        dst_dir_4ch = fullfile(dst_dir, view_str);
        mkdir(dst_dir_4ch);

        gt_4ch_norm = gt_4ch;

        RO = size(gt_4ch,1);
        E1 = size(gt_4ch,2);
        ps = max(gt_4ch_ps);

        new_RO = round(RO*ps/dst_pixel_spacing(1));
        new_E1 = round(E1*ps/dst_pixel_spacing(2));

        data = Matlab_gt_resize_2D_image(double(gt_4ch), new_RO, new_E1, 5);

        if(plot_flag)
            figure; imagescn(data, [], [], [8], 3);
        end

        data_dst = data;
        RO = size(data_dst,1);
        E1 = size(data_dst,2);
        PHS = size(data_dst,3);

        data = zpad(data,max(RO, NN_RO), max(E1, NN_E1), PHS);
        data_dst = zpad(data_dst,max(RO, NN_RO), max(E1, NN_E1), PHS);

        RO = size(data_dst,1);
        E1 = size(data_dst,2);
        PHS = size(data_dst,3);

        [mask, ro_s, ro_e, e1_s, e1_e, h_mask] = detect_heart_region_cine_roi(data_dst, 0.9, NN_RO, NN_E1, use_moco, regularization_hilbert_strength);
        set(h_mask, 'visible', 'off');

        Cine_resized_training = data; % data(ro_s:ro_e, e1_s:e1_e, :, :);
        h2 = figure('visible', visible_status); imagescn(Cine_resized_training, [], [], [8], 3);
        set(h2, 'visible', 'off')

        header = CreateGtImageHeader(Cine_resized_training);
        Matlab_gt_write_analyze(single(Cine_resized_training), header, fullfile(dst_dir_4ch, 'Cine_resized_training'));
        writeNPY(single(gt_4ch), fullfile(dst_dir_4ch, 'data_acq.npy'));
        writeNPY(single(Cine_resized_training), fullfile(dst_dir_4ch, 'Cine_resized_training.npy'));

        writeNPY(single(data), fullfile(dst_dir_4ch, 'data.npy'));
        writeNPY(single([ro_s, ro_e, e1_s, e1_e]), fullfile(dst_dir_4ch, 'roi.npy'));

        saveas(h_mask, fullfile(dst_dir_4ch, 'roi'), 'jpg');
        saveas(h2, fullfile(dst_dir_4ch, 'Cine_resized_training'), 'jpg');

        pixel_spacing = gt_4ch_ps;
        save(fullfile(dst_dir_4ch, 'record'), 'dst_pixel_spacing', 'pixel_spacing');
        save(fullfile(dst_dir_4ch, 'record_header'), 'dst_pixel_spacing', 'pixel_spacing');

        for p=1:numel(phases)
            if(phases(p)>PHS)
                continue;
            end
            im = data_dst(:,:,phases(p));

%             ED_im = 255*im/max(im(:));
%             ED_im = normalizeWindowSetting(ED_im, 3.0*mean(ED_im(:)), 4.0*mean(ED_im(:)));
%             ED_im = uint8(ED_im);
            
            L = im/max(im(:));
            L2 = adapthisteq(L,'NumTiles',[8 8],'ClipLimit',0.005);
            ED_im = uint8(255*L2);
            
            imwrite(ED_im, fullfile(dst_dir_4ch, [case_prefix '_' view_str '_' num2str(phases(p)) '.jpg']), 'jpg');
            imwrite(ED_im, fullfile(pic_dir, [case_prefix '_' view_str '_' num2str(phases(p)) '.jpg']), 'jpg');

            writeNPY(single(im), fullfile(pic_dir_numpy, [case_prefix '_' view_str '_' num2str(phases(p)) '.npy']));
        end
    catch
        disp(['Exceptions happened for ' case_prefix]);
    end
end
