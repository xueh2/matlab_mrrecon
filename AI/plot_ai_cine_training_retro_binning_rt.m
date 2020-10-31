function plot_ai_cine_training_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, linewidth, ED, ES, check_processed )
% plot_ai_cine_training_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, linewidth, ED, ES)

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
        if(~check_processed | ~exist(fullfile(contourDir, 'NN_contours.mat')))
            
            record = load(fullfile(dst_dir_sax, 'record'));
            [endo, epi, Cine_resized_training] = cine_load_NN_res_on_images(dst_dir_sax, contourDir, load_endo, load_epi);
            
            PHS = size(record.gt_sax_norm, 4);
            SLC = size(record.gt_sax_norm, 3);
            
            slc_loc = zeros(size(record.gt_sax_norm, 3),1);
            for slc=1:size(record.gt_sax_norm, 3)
                slc_loc(slc) = record.gt_h_sax(slc).slice_location;
            end
            
            [sv, ind_slc] = sort(slc_loc);
            
            start_slc = floor(SLC/2)-1;
            end_slc = ceil(SLC/2)+1;
            
            phs_lv_area = zeros(PHS, 1);
            
            for tt=start_slc:end_slc
                slc = ind_slc(tt);
                good_slc = 1;
                for phs=1:PHS
                    endo_seg = endo.seg(:,:,phs,slc);
                    if(sum(endo_seg(:))<=10)
                        good_slc = 0;
                        break;
                    end
                end
                
                if(good_slc)
                    for phs=1:PHS
                        endo_seg = endo.seg(:,:,phs,slc);
                        phs_lv_area(phs) =  phs_lv_area(phs) + sum(endo_seg(:));
                    end
                end
            end
            
            [min_lv_area, ES] = min(phs_lv_area);
            
            h = figure;
            hold on
            plot(1:PHS, phs_lv_area);
            hold off
            box on
            title('lv area');                        

            % plot contours
            % [h_slc, h_ED, h_ES, h_slc_obj] = cine_plot_contours_on_images(Cine_resized_training, endo.C, epi.C, linewidth, ED, ES);
            [h_slc, h_ED, h_ES, h_slc_obj] = cine_plot_contours_on_images(Cine_resized_training, endo.contour_C, epi.contour_C, linewidth, ED, ES);

            for slc=1:numel(h_slc)
                saveas(h_slc(slc), fullfile(contourDir, ['SLC_' num2str(slc)]), 'fig');
                saveas(h_slc_obj(slc), fullfile(contourDir, ['SLC_obj_' num2str(slc)]), 'fig');
            end

            saveas(h_ED, fullfile(contourDir, ['ED']), 'fig');
            saveas(h_ES, fullfile(contourDir, ['ES']), 'fig');

            save(fullfile(contourDir, 'NN_contours'), 'endo', 'epi', 'ES', 'phs_lv_area');
            saveas(h, fullfile(contourDir, 'lv_area'), 'jpg');
            
            cine_add_endo_epi_contours(Cine_resized_training, endo.contour_C, epi.contour_C, ED, contourDir);
            cine_add_endo_epi_contours(Cine_resized_training, endo.contour_C, epi.contour_C, ES, contourDir);
        end
    end
end
