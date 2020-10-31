function prepare_ai_cine_training_retro_binning_rt_making_jpg_pics(resDir, aiDir, pt_ids, files_record_picked, check_processed)
% prepare_ai_cine_training_retro_binning_rt_making_jpg_pics(resDir, aiDir, pt_ids, files_record_picked, check_processed)

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

    if(numel(case_4ch)==0 || numel(case_2ch)==0 || numel(case_sax)==0)
        continue;
    end
        
    case_prefix = [case_4ch.study_dates(end,:) '_' pt_id];
    dst_dir = fullfile(aiDir, case_4ch.study_dates(end,:), pt_id);
    
    if(check_processed)
        dst_dir_sax = fullfile(dst_dir, 'sax');
        if(numel(case_3ch)==0)
            if(exist(fullfile(dst_dir_sax, 'record_header.mat')) ...
                & exist(fullfile(pic_dir, 'ch2', [case_prefix '_ch2_' 'ES.jpg'])) ...
                & exist(fullfile(pic_dir, 'ch4', [case_prefix '_ch4_' 'ES.jpg'])) ...
                & exist(fullfile(pic_dir, 'sax', [case_prefix '_sax_' 'ES_6.jpg'])) ...
                & exist(fullfile(pic_dir, 'sax_original_numpy', [case_prefix '_sax_' 'ED_6.npy'])) )
            
                disp(['already processed ' pt_id]);

                continue;
            end
        else
            if(exist(fullfile(dst_dir_sax, 'record_header.mat')) ...
                    & exist(fullfile(pic_dir, 'ch2', [case_prefix '_ch2_' 'ES.jpg'])) ...
                    & exist(fullfile(pic_dir, 'ch3', [case_prefix '_ch3_' 'ES.jpg'])) ...
                    & exist(fullfile(pic_dir, 'ch4', [case_prefix '_ch4_' 'ES.jpg'])) ...
                    & exist(fullfile(pic_dir, 'sax', [case_prefix '_sax_' 'ES_6.jpg'])) ...
                    & exist(fullfile(pic_dir, 'sax_original_numpy', [case_prefix '_sax_' 'ED_6.npy'])) )

                disp(['already processed ' pt_id]);

                continue;
            end
        end
    end
    
    SLC = case_sax.headers{end}.encoding.encodingLimits.slice.maximum +1;
    sax_ind = size(case_sax, 1);
    if(SLC<6)
        SLC = case_sax.headers{1}.encoding.encodingLimits.slice.maximum +1;            
        sax_ind = 1;
    end

    % -----------------------------

    case_4ch_dir = fullfile(resDir, case_4ch.study_dates(end,:), case_4ch.file_names{end})

    % -----------------------------

    case_3ch_dir = fullfile(resDir, case_3ch.study_dates(end,:), case_3ch.file_names{end})

    % -----------------------------

    case_2ch_dir = fullfile(resDir, case_2ch.study_dates(end,:), case_2ch.file_names{end})       
    case_sax_dir = fullfile(resDir, case_sax.study_dates(sax_ind,:), case_sax.file_names{sax_ind})       
    
    mkdir(dst_dir)

    % 4ch
    dst_dir_4ch = fullfile(dst_dir, 'ch4');
       
    try
        data = readNPY(fullfile(dst_dir_4ch, 'data.npy'));
        data_dst = readNPY(fullfile(dst_dir_4ch, 'data_norm.npy'));
    catch
        continue;
    end
    
    ED_phs = data_dst(:,:,ED);
    ES_phs = data_dst(:,:,ES);
    
    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_ch4_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_ch4_' 'ES.jpg']), 'jpg');
    
    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(pic_dir, 'ch4', [case_prefix '_ch4_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(pic_dir, 'ch4', [case_prefix '_ch4_' 'ES.jpg']), 'jpg');
    
    writeNPY(single(ED_phs), fullfile(pic_dir, 'ch4_numpy', [case_prefix '_ch4_' 'ED.npy']));
    writeNPY(single(ES_phs), fullfile(pic_dir, 'ch4_numpy', [case_prefix '_ch4_' 'ES.npy']));
    
    ED_phs = data(:,:,ED);
    ES_phs = data(:,:,ES);
    
    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_ch4_original_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(dst_dir_4ch, [case_prefix '_ch4_original_' 'ES.jpg']), 'jpg');
    
    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(pic_dir, 'ch4_original', [case_prefix '_ch4_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(pic_dir, 'ch4_original', [case_prefix '_ch4_' 'ES.jpg']), 'jpg');
    
    writeNPY(single(ED_phs), fullfile(pic_dir, 'ch4_original_numpy', [case_prefix '_ch4_' 'ED.npy']));
    writeNPY(single(ES_phs), fullfile(pic_dir, 'ch4_original_numpy', [case_prefix '_ch4_' 'ES.npy']));
    
    % 2ch
    dst_dir_2ch = fullfile(dst_dir, 'ch2');

    data = readNPY(fullfile(dst_dir_2ch, 'data.npy'));
    data_dst = readNPY(fullfile(dst_dir_2ch, 'data_norm.npy'));
    
    ED_phs = data_dst(:,:,ED);
    ES_phs = data_dst(:,:,ES);
    
    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(dst_dir_2ch, [case_prefix '_ch2_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(dst_dir_2ch, [case_prefix '_ch2_' 'ES.jpg']), 'jpg');
    
    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(pic_dir, 'ch2', [case_prefix '_ch2_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(pic_dir, 'ch2', [case_prefix '_ch2_' 'ES.jpg']), 'jpg');
    
    writeNPY(single(ED_phs), fullfile(pic_dir, 'ch2_numpy', [case_prefix '_ch2_' 'ED.npy']));
    writeNPY(single(ES_phs), fullfile(pic_dir, 'ch2_numpy', [case_prefix '_ch2_' 'ES.npy']));
    
    ED_phs = data(:,:,ED);
    ES_phs = data(:,:,ES);
    
    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(dst_dir_2ch, [case_prefix '_ch2_original_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(dst_dir_2ch, [case_prefix '_ch2_original_' 'ES.jpg']), 'jpg');
    
    imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(pic_dir, 'ch2_original', [case_prefix '_ch2_' 'ED.jpg']), 'jpg');
    imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(pic_dir, 'ch2_original', [case_prefix '_ch2_' 'ES.jpg']), 'jpg');
    
    writeNPY(single(ED_phs), fullfile(pic_dir, 'ch2_original_numpy', [case_prefix '_ch2_' 'ED.npy']));
    writeNPY(single(ES_phs), fullfile(pic_dir, 'ch2_original_numpy', [case_prefix '_ch2_' 'ES.npy']));
    
    % 3ch
    if(numel(case_3ch)>0)
        dst_dir_3ch = fullfile(dst_dir, 'ch3');
            
        data = readNPY(fullfile(dst_dir_3ch, 'data.npy'));
        data_dst = readNPY(fullfile(dst_dir_3ch, 'data_norm.npy'));

        ED_phs = data_dst(:,:,ED);
        ES_phs = data_dst(:,:,ES);

        imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(dst_dir_3ch, [case_prefix '_ch3_' 'ED.jpg']), 'jpg');
        imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(dst_dir_3ch, [case_prefix '_ch3_' 'ES.jpg']), 'jpg');

        imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(pic_dir, 'ch3', [case_prefix '_ch3_' 'ED.jpg']), 'jpg');
        imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(pic_dir, 'ch3', [case_prefix '_ch3_' 'ES.jpg']), 'jpg');

        writeNPY(single(ED_phs), fullfile(pic_dir, 'ch3_numpy', [case_prefix '_ch3_' 'ED.npy']));
        writeNPY(single(ES_phs), fullfile(pic_dir, 'ch3_numpy', [case_prefix '_ch3_' 'ES.npy']));

        ED_phs = data(:,:,ED);
        ES_phs = data(:,:,ES);

        imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(dst_dir_3ch, [case_prefix '_ch3_original_' 'ED.jpg']), 'jpg');
        imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(dst_dir_3ch, [case_prefix '_ch3_original_' 'ES.jpg']), 'jpg');

        imwrite(uint8(255*ED_phs/max(ED_phs(:))), fullfile(pic_dir, 'ch3_original', [case_prefix '_ch3_' 'ED.jpg']), 'jpg');
        imwrite(uint8(255*ES_phs/max(ES_phs(:))), fullfile(pic_dir, 'ch3_original', [case_prefix '_ch3_' 'ES.jpg']), 'jpg');

        writeNPY(single(ED_phs), fullfile(pic_dir, 'ch3_original_numpy', [case_prefix '_ch3_' 'ED.npy']));
        writeNPY(single(ES_phs), fullfile(pic_dir, 'ch3_original_numpy', [case_prefix '_ch3_' 'ES.npy']));
    end
    
    %sax
    dst_dir_sax = fullfile(dst_dir, 'sax');
            
    data = readNPY(fullfile(dst_dir_sax, 'data.npy'));
    data_dst = readNPY(fullfile(dst_dir_sax, 'data_norm.npy'));

    ED_phs = squeeze(data_dst(:,:,ED,:));
    ES_phs = squeeze(data_dst(:,:,ES,:));
    
    for slc=1:SLC
        im = ED_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(dst_dir_sax, [case_prefix '_sax_' 'ED_' num2str(slc) '.jpg']), 'jpg');
        im = ES_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(dst_dir_sax, [case_prefix '_sax_' 'ES_' num2str(slc) '.jpg']), 'jpg');
    end
    
    for slc=1:SLC
        im = ED_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(pic_dir, 'sax', [case_prefix '_sax_' 'ED_' num2str(slc) '.jpg']), 'jpg');
        im = ES_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(pic_dir, 'sax', [case_prefix '_sax_' 'ES_' num2str(slc) '.jpg']), 'jpg');
        
        writeNPY(single(ED_phs(:,:,slc)), fullfile(pic_dir, 'sax_numpy', [case_prefix '_sax_' 'ED_' num2str(slc) '.npy']));
        writeNPY(single(ES_phs(:,:,slc)), fullfile(pic_dir, 'sax_numpy', [case_prefix '_sax_' 'ES_' num2str(slc) '.npy']));
    end
    
    ED_phs = squeeze(data(:,:,ED,:));
    ES_phs = squeeze(data(:,:,ES,:));
    
    for slc=1:SLC
        im = ED_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(dst_dir_sax, [case_prefix '_sax_original_' 'ED_' num2str(slc) '.jpg']), 'jpg');
        im = ES_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(dst_dir_sax, [case_prefix '_sax_original_' 'ES_' num2str(slc) '.jpg']), 'jpg');
    end
    
    for slc=1:SLC
        im = ED_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(pic_dir, 'sax_original', [case_prefix '_sax_' 'ED_' num2str(slc) '.jpg']), 'jpg');
        im = ES_phs(:,:,slc);
        imwrite(uint8(255*im/max(im(:))), fullfile(pic_dir, 'sax_original', [case_prefix '_sax_' 'ES_' num2str(slc) '.jpg']), 'jpg');
        
        writeNPY(single(ED_phs(:,:,slc)), fullfile(pic_dir, 'sax_original_numpy', [case_prefix '_sax_' 'ED_' num2str(slc) '.npy']));
        writeNPY(single(ES_phs(:,:,slc)), fullfile(pic_dir, 'sax_original_numpy', [case_prefix '_sax_' 'ES_' num2str(slc) '.npy']));
    end
    
    closeall
    closeall    
end
