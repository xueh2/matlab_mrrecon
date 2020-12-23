function collect_ai_perf_pd(dataDir, resDir, aiDir, pt_ids, files_record_picked, suffix, mon, checkprocessed)
% collect_ai_perf_pd(dataDir, resDir, aiDir, pt_ids, files_record_picked, suffix, mon, checkprocessed)

ai_res_dir = fullfile(aiDir, 'ai_results_for_training');
mkdir(ai_res_dir);

stress_pic_dir = fullfile(ai_res_dir, ['stress_month_' num2str(mon)]);
mkdir(stress_pic_dir);
stress_fid_all = fopen(fullfile(ai_res_dir, ['stress_ai_res_month_' num2str(mon) '.csv']), 'w+');
fprintf(stress_fid_all, '%s\n', 'filename,file_size,file_attributes,region_count,region_id,region_shape_attributes,region_attributes');

rest_pic_dir = fullfile(ai_res_dir, ['rest_month_' num2str(mon)]);
mkdir(rest_pic_dir);
rest_fid_all = fopen(fullfile(ai_res_dir, ['rest_ai_res_month_' num2str(mon) '.csv']), 'w+');
fprintf(rest_fid_all, '%s\n', 'filename,file_size,file_attributes,region_count,region_id,region_shape_attributes,region_attributes');

stress_moco_pic_dir = fullfile(ai_res_dir, ['stress_moco_month_' num2str(mon)]);
mkdir(stress_moco_pic_dir);
stress_moco_fid_all = fopen(fullfile(ai_res_dir, ['stress_moco_ai_res_month_' num2str(mon) '.csv']), 'w+');
fprintf(stress_moco_fid_all, '%s\n', 'filename,file_size,file_attributes,region_count,region_id,region_shape_attributes,region_attributes');

rest_moco_pic_dir = fullfile(ai_res_dir, ['rest_moco_month_' num2str(mon)]);
mkdir(rest_moco_pic_dir);
rest_moco_fid_all = fopen(fullfile(ai_res_dir, ['rest_moco_ai_res_month_' num2str(mon) '.csv']), 'w+');
fprintf(rest_moco_fid_all, '%s\n', 'filename,file_size,file_attributes,region_count,region_id,region_shape_attributes,region_attributes');

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    if(isnumeric(pt_id))
        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' num2str(pt_id)]);
    else
        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    end
    
    case_stress = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'stress', 1);
    case_rest = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'rest', 1);
       
    if(numel(case_stress)>0)
        study_date = case_stress.study_dates(end,:);
        pt_id = case_stress.patientIDs{end,:};
    end
    if(numel(case_rest)>0)
        study_date = case_rest.study_dates(end,:);
        pt_id = case_rest.patientIDs{end,:};
    end
    disp([pt_id ' - ' study_date]);
           
    dst_dir = fullfile(aiDir, study_date, pt_id);
    if(~exist(dst_dir))
        dst_dir = fullfile(aiDir, pt_id);
    end
    
    % stress
    dst_dir_stress = fullfile(dst_dir, 'stress');
               
    % rest
    dst_dir_rest = fullfile(dst_dir, 'rest');
           
    contourDir = fullfile(dst_dir, 'res_ai');
          
    if(exist(fullfile(contourDir, ['stress_pd_seg_C' suffix '.npy'])))
        
        case_name = case_stress.file_names{1};
        
        pd = readNPY(fullfile(dst_dir_stress, 'pd_for_seg.npy'));
        masks = readNPY(fullfile(contourDir, ['stress_pd_seg_masks' suffix '.npy']));
        contours = readNPY(fullfile(contourDir, ['stress_pd_seg_C' suffix '.npy']));

        process_one_case(pd, contours, contourDir, stress_pic_dir, case_name, 'stress_pd_seg_C.csv', stress_fid_all);
    end        
    
    if(exist(fullfile(contourDir, ['rest_pd_seg_C' suffix '.npy'])))
        
        case_name = case_rest.file_names{1};
        
        pd = readNPY(fullfile(dst_dir_rest, 'pd_for_seg.npy'));
        masks = readNPY(fullfile(contourDir, ['rest_pd_seg_masks' suffix '.npy']));
        contours = readNPY(fullfile(contourDir, ['rest_pd_seg_C' suffix '.npy']));

        process_one_case(pd, contours, contourDir, rest_pic_dir, case_name, 'rest_pd_seg_C.csv', rest_fid_all);
    end 
    
    if(exist(fullfile(contourDir, ['stress_pd_moco_seg_C' suffix '.npy'])))
        
        case_name = case_stress.file_names{1};
        
        pd = readNPY(fullfile(dst_dir_stress, 'pd_moco_for_seg.npy'));
        masks = readNPY(fullfile(contourDir, ['stress_pd_moco_seg_masks' suffix '.npy']));
        contours = readNPY(fullfile(contourDir, ['stress_pd_moco_seg_C' suffix '.npy']));

        process_one_case(pd, contours, contourDir, stress_moco_pic_dir, case_name, 'stress_pd_moco_seg_C.csv', stress_moco_fid_all);
    end
    if(exist(fullfile(contourDir, ['rest_pd_moco_seg_C' suffix '.npy'])))
        
        case_name = case_rest.file_names{1};
        
        pd = readNPY(fullfile(dst_dir_rest, 'pd_moco_for_seg.npy'));
        masks = readNPY(fullfile(contourDir, ['rest_pd_moco_seg_masks' suffix '.npy']));
        contours = readNPY(fullfile(contourDir, ['rest_pd_moco_seg_C' suffix '.npy']));

        process_one_case(pd, contours, contourDir, rest_moco_pic_dir, case_name, 'rest_pd_moco_seg_C.csv', rest_moco_fid_all);
    end
end
fclose(stress_fid_all);
fclose(rest_fid_all);
fclose(stress_moco_fid_all);
fclose(rest_moco_fid_all);

end

function process_one_case(pd, contours, contourDir, pic_dir, case_name, csv_name, fid_all)
    pd = squeeze(pd);
    SLC = size(pd, 3);

    C = [];
    for slc=1:SLC
        pd_im = pd(:,:,slc);
        mean_v = mean(pd_im(:));
        im = uint8(255 * pd_im/(3 * mean_v));
        imwrite(im, fullfile(pic_dir, [case_name '_slc_' num2str(slc) '.jpg']), 'jpg');            
        C{slc} = contours(:,[2 1],slc) + 1; %mask2contour(masks(:,:,slc), 1, 32, 32);
    end

    fid = fopen(fullfile(contourDir, csv_name), 'w+');
    for slc=1:SLC
        info = dir(fullfile(pic_dir, [case_name '_slc_' num2str(slc) '.jpg']));
        curr_C = C{slc};            
        csv_str = generate_VIA_polygon_str(curr_C, [case_name '_slc_' num2str(slc) '.jpg'], info);
        fprintf(fid, '%s\n', csv_str);
        fprintf(fid_all, '%s\n', csv_str);
    end
    fclose(fid);
end
