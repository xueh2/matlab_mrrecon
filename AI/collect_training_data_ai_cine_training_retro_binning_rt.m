function collect_training_data_ai_cine_training_retro_binning_rt(resDir, aiDir, trainDir_used, pt_ids, files_record_picked)
% collect_training_data_ai_cine_training_retro_binning_rt(resDir, aiDir, trainDir_used, pt_ids, files_record_picked)

load_endo = 1;
load_epi = 1;

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    
    case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '4ch');
    case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '2ch');
    case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'sa');

    if(numel(case_4ch)==0 || numel(case_2ch)==0 || numel(case_sax)==0)
        continue;
    end
        
    dst_dir = fullfile(aiDir, case_4ch.study_dates(end,:), pt_id);
    
    try
        SLC = case_sax.headers{end}.encoding.encodingLimits.slice.maximum +1;
        sax_ind = size(case_sax, 1);
        if(SLC<6)
            SLC = case_sax.headers{1}.encoding.encodingLimits.slice.maximum +1;            
            sax_ind = 1;
        end

        case_4ch_dir = fullfile(resDir, case_4ch.study_dates(end,:), case_4ch.file_names{end})       
        case_2ch_dir = fullfile(resDir, case_2ch.study_dates(end,:), case_2ch.file_names{end})       
        case_sax_dir = fullfile(resDir, case_sax.study_dates(end,:), case_sax.file_names{sax_ind})       
    catch
        continue;
    end
    
    % 4ch
    dst_dir_4ch = fullfile(dst_dir, 'ch4');
               
    % 2ch
    dst_dir_2ch = fullfile(dst_dir, 'ch2');
           
    %sax
    dst_dir_sax = fullfile(dst_dir, 'sax');
       
    contourDir = fullfile(dst_dir, 'sax_ai');
   
    if(exist(fullfile(contourDir, 'Cine_resized_training_5.hdr')))    
        if(exist(fullfile(contourDir, 'NN_contours.mat')))

            roi_files = findFILE(contourDir, '*new.mat');
            
            if(numel(roi_files)==0)
                continue;
            end
            
            Cine_resized_training = readNPY(fullfile(dst_dir_sax, 'Cine_resized_training_norm.npy'));
            im_roi = readNPY(fullfile(dst_dir_sax, 'roi.npy'));
            data = readNPY(fullfile(dst_dir_sax, 'data.npy'));
            data_norm = readNPY(fullfile(dst_dir_sax, 'data_norm.npy'));
         
            RO = size(data, 1);
            E1 = size(data, 2);
            PHS = size(data, 3);
            SLC = size(data, 4);
            
            for kk=1:numel(roi_files)
            
                roi_name = roi_files{kk};
                disp([case_sax.study_dates(end,:) ' - ' case_sax.file_names{sax_ind} ' - ' roi_name]);
                
                ind = strfind(roi_name, 'roi_phs');
                ind2 = strfind(roi_name, '_new');                
                phs = str2double(roi_name(ind+7:ind2-1));
                
                roi = load(roi_name);
                
                ED = phs;
                plotFlag = 1;
                [S_PHS, h] = cine_generate_seg_from_manual_roi(Cine_resized_training, roi.ROI_info_table, ED, plotFlag);

                train_data_dir = fullfile(trainDir_used, case_4ch.study_dates(end,:), pt_id);
                mkdir(train_data_dir);
                
                endo_mask_full_fov = zeros(RO, E1, SLC);
                epi_mask_full_fov = zeros(RO, E1, SLC);
                rvi_mask_full_fov = zeros(RO, E1, SLC);
                
                endo_mask_full_fov(im_roi(1):im_roi(2), im_roi(3):im_roi(4), :) = S_PHS.endo_mask;
                epi_mask_full_fov(im_roi(1):im_roi(2), im_roi(3):im_roi(4), :) = S_PHS.epi_mask;
                rvi_mask_full_fov(im_roi(1):im_roi(2), im_roi(3):im_roi(4), :) = S_PHS.rvi_mask;
                
                SLC = size(S_PHS.im, 3);
                
                for slc=1:SLC
                    train_slc_dir = fullfile(train_data_dir, ['slc_' num2str(slc) '_phs_' num2str(phs)]);
                    mkdir(train_slc_dir);
                    
                    im = data(:,:,phs,slc);
                    im_norm = data_norm(:,:,phs,slc);
                    endo = endo_mask_full_fov(:,:,slc);
                    epi = epi_mask_full_fov(:,:,slc);
                    rvi = rvi_mask_full_fov(:,:,slc);
                                                            
                    writeNPY(single(im), fullfile(train_slc_dir, 'im_train.npy'));
                    writeNPY(single(im_norm), fullfile(train_slc_dir, 'im_norm_train.npy'));
                    writeNPY(single(endo), fullfile(train_slc_dir, 'endo_train.npy'));
                    writeNPY(single(epi), fullfile(train_slc_dir, 'epi_train.npy'));
                    writeNPY(single(rvi), fullfile(train_slc_dir, 'rvi_train.npy'));                                        
                    writeNPY(single(im_roi), fullfile(train_slc_dir, 'roi.npy'));     
                    
                    header = CreateGtImageHeader(im);
                    Matlab_gt_write_analyze(single(im), header, fullfile(train_slc_dir, 'im'));
                    Matlab_gt_write_analyze(single(im_norm), header, fullfile(train_slc_dir, 'im_norm'));
                                        
                    h2 = figure;
                    imagescn(cat(3, im, im_norm, endo, epi, rvi), [], [1 5], [12]);
                    saveas(h2, fullfile(train_slc_dir, 'im_train'), 'jpg');
                end
                
                saveas(h, fullfile(train_data_dir, ['cine_train_phs' num2str(phs)]), 'jpg');
                
                closeall
                closeall
            end
        end
    end
end
