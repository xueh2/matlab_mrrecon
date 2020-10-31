function plot_dicom_training_data_ai_cine_training_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, series_num, check_processed)
% plot_dicom_training_data_ai_cine_training_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, series_num, check_processed)
% expect input as column-wise GT saved images (output_image_in_column_storage=true)

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
    if(exist(dst_dir_4ch)==0)
        continue;
    end
    
    % 2ch
    dst_dir_2ch = fullfile(dst_dir, 'ch2');
    if(exist(dst_dir_2ch)==0)
        continue;
    end
    
    %sax
    dst_dir_sax = fullfile(dst_dir, 'sax');
    if(exist(dst_dir_sax)==0)
        continue;
    end
    
    contourDir = fullfile(dst_dir, 'sax_ai');

    if(check_processed & exist(fullfile(dst_dir, 'sax_1mm_pixel_spacing.npy')))
        continue;
    end
    
    try
        gt_4ch_norm = readNPY(fullfile(dst_dir_4ch, 'data_acq_norm.npy'));
        gt_2ch_norm = readNPY(fullfile(dst_dir_2ch, 'data_acq_norm.npy'));
        gt_sax_norm = readNPY(fullfile(dst_dir_sax, 'data_acq_norm.npy'));
        gt_4ch = gt_4ch_norm;
        gt_2ch = gt_2ch_norm;
        gt_sax = gt_sax_norm;
        
        gt_4ch_norm_1mm = readNPY(fullfile(dst_dir_4ch, 'data_norm.npy'));
        gt_2ch_norm_1mm = readNPY(fullfile(dst_dir_2ch, 'data_norm.npy'));
        gt_sax_norm_1mm = readNPY(fullfile(dst_dir_sax, 'data_norm.npy'));
        
        gt_4ch = permute(gt_4ch, [2 1 3]);
        gt_2ch = permute(gt_2ch, [2 1 3]);
        gt_sax_norm = permute(gt_sax_norm, [2 1 4 3]);
        
        gt_4ch_norm_1mm = permute(gt_4ch_norm_1mm, [2 1 3]);
        gt_2ch_norm_1mm = permute(gt_2ch_norm_1mm, [2 1 3]);
        gt_sax_norm_1mm = permute(gt_sax_norm_1mm, [2 1 4 3]);
    catch
        continue;
    end
    
    record = load(fullfile(dst_dir_4ch, 'record_header'));
    gt_h_4ch = record.gt_h_4ch;
    
    record = load(fullfile(dst_dir_2ch, 'record_header'));
    gt_h_2ch = record.gt_h_2ch;
    
    record = load(fullfile(dst_dir_sax, 'record_header'));
    gt_h_sax = record.gt_h_sax;
    gt_h_sax = squeeze(gt_h_sax);
    
    % make image in cmr view for dicom convention
    rotate90_4ch= -1;
    flip_4ch = 2;        
    [I_4ch, gt_4ch_h_cmr] = apply_normorientation(gt_4ch(:,:,1), rotate90_4ch, flip_4ch, gt_h_4ch(1));
    h_gt_4ch_cmr = ComputeDicomCoordFromGtOffline(gt_4ch_h_cmr(1).PatientPosition, gt_4ch_h_cmr(1).read_dir, gt_4ch_h_cmr(1).phase_dir, double(gt_4ch_h_cmr(1).FOV), double([size(I_4ch, 1), size(I_4ch, 2), 1]), size(I_4ch, 2), size(I_4ch, 1));    

    rotate90_2ch= -1;
    flip_2ch = 2;        
    [I_2ch, gt_2ch_h_cmr] = apply_normorientation(gt_2ch(:,:,1), rotate90_2ch, flip_2ch, gt_h_2ch(1));
    gt_2ch_h = gt_h_2ch(1, 1);
    h_gt_2ch_cmr = ComputeDicomCoordFromGtOffline(gt_2ch_h_cmr(1).PatientPosition, gt_2ch_h_cmr(1).read_dir, gt_2ch_h_cmr(1).phase_dir, double(gt_2ch_h_cmr(1).FOV), double([size(I_2ch, 1), size(I_2ch, 2), 1]), size(I_2ch, 2), size(I_2ch, 1));

    clear I_sax gt_sax_h_cmr  h_gt_sax_cmr
    SLC = size(gt_sax_norm_1mm, 3);
    
    for slc=1:SLC    
        rotate90_sax= -1;
        flip_sax = 2    ;    
        [I_sax(:,:,slc), gt_sax_h_cmr(slc)] = apply_normorientation(gt_sax_norm(:,:,slc,1), rotate90_sax, flip_sax, gt_h_sax(slc));
        h_gt_sax_cmr(slc) = ComputeDicomCoordFromGtOffline(gt_sax_h_cmr(slc).PatientPosition, gt_sax_h_cmr(slc).read_dir, gt_sax_h_cmr(slc).phase_dir, double(gt_sax_h_cmr(slc).FOV), double([size(I_sax, 1), size(I_sax, 2), 1]), size(I_sax, 2), size(I_sax, 1));
    end
    
    % -------------------------------------
    h1 = figure;
    hold on
    ah = axes;
    set(gcf, 'Renderer', 'zbuffer');
    plot2DDicomSlice(I_4ch, h_gt_4ch_cmr, ah, -1);
    plot2DDicomSlice(I_2ch, h_gt_2ch_cmr, ah, -1);
    if(SLC>5)
        plot2DDicomSlice(I_sax(:,:,5), h_gt_sax_cmr(5), ah, -1);
    else
        plot2DDicomSlice(I_sax(:,:,SLC), h_gt_sax_cmr(SLC), ah, -1);
    end
    I_4ch = zpad(I_4ch, [512, 512]);
    I_2ch = zpad(I_2ch, [512, 512]);
    I_sax = zpad(I_sax, 512, 512, SLC);
    
    if(SLC>5)
        h2 = figure; imagescn(cat(3, I_4ch, I_2ch, I_sax(:,:,5)), [], [1 3], 12);
    else
        h2 = figure; imagescn(cat(3, I_4ch, I_2ch, I_sax(:,:,SLC)), [], [1 3], 12);
    end
    
    saveas(h1, fullfile(dst_dir, 'multi_slice_view'), 'jpg');
    saveas(h2, fullfile(dst_dir, 'ch4_ch2_sax_in_one'), 'jpg');
    
    % -------------------------------------
    I_4ch_all = permute(gt_4ch, [2 1 3]);
    I_2ch_all = permute(gt_2ch, [2 1 3]);
    I_sax_all = permute(gt_sax_norm, [2 1 4 3]);
    
    I_4ch_1mm_all = permute(gt_4ch_norm_1mm, [2 1 3]);
    I_2ch_1mm_all = permute(gt_2ch_norm_1mm, [2 1 3]);
    I_sax_1mm_all = permute(gt_sax_norm_1mm, [2 1 4 3]);
    
    writeNPY(single(I_4ch_all), fullfile(dst_dir, 'I_4ch.npy'));
    writeNPY(single(I_2ch_all), fullfile(dst_dir, 'I_2ch.npy'));
    writeNPY(single(I_sax_all), fullfile(dst_dir, 'I_sax.npy'));
    
    writeNPY(single(I_4ch_1mm_all), fullfile(dst_dir, 'I_4ch_1mm.npy'));
    writeNPY(single(I_2ch_1mm_all), fullfile(dst_dir, 'I_2ch_1mm.npy'));
    writeNPY(single(I_sax_1mm_all), fullfile(dst_dir, 'I_sax_1mm.npy'));
    
    Matlab_gt_write_analyze(single(I_4ch_all), CreateGtImageHeader(I_4ch_all), fullfile(dst_dir, 'I_4ch'));
    Matlab_gt_write_analyze(single(I_2ch_all), CreateGtImageHeader(I_2ch_all), fullfile(dst_dir, 'I_2ch'));
    Matlab_gt_write_analyze(single(I_sax_all), CreateGtImageHeader(I_sax_all), fullfile(dst_dir, 'I_sax'));
    
    Matlab_gt_write_analyze(single(I_4ch_1mm_all), CreateGtImageHeader(I_4ch_all), fullfile(dst_dir, 'I_4ch_1mm'));
    Matlab_gt_write_analyze(single(I_2ch_1mm_all), CreateGtImageHeader(I_2ch_all), fullfile(dst_dir, 'I_2ch_1mm'));
    Matlab_gt_write_analyze(single(I_sax_1mm_all), CreateGtImageHeader(I_sax_all), fullfile(dst_dir, 'I_sax_1mm'));
    
    writeNPY(single(h_gt_4ch_cmr.ImagePositionPatient), fullfile(dst_dir, 'ch4_ImagePosition.npy'));
    writeNPY(single(h_gt_4ch_cmr.ImageOrientationPatient), fullfile(dst_dir, 'ch4_ImageOrientation.npy'));
    writeNPY(single([h_gt_4ch_cmr.PixelSpacing h_gt_4ch_cmr.SliceThickness]), fullfile(dst_dir, 'ch4_pixel_spacing.npy'));
    writeNPY(single([1 1 h_gt_4ch_cmr.SliceThickness]), fullfile(dst_dir, 'ch4_1mm_pixel_spacing.npy'));

    Matlab_gt_write_analyze(single(h_gt_4ch_cmr.ImagePositionPatient), CreateGtImageHeader(h_gt_4ch_cmr.ImagePositionPatient), fullfile(dst_dir, 'ch4_ImagePosition'));
    Matlab_gt_write_analyze(single(h_gt_4ch_cmr.ImageOrientationPatient), CreateGtImageHeader(h_gt_4ch_cmr.ImageOrientationPatient), fullfile(dst_dir, 'ch4_ImageOrientation'));
    Matlab_gt_write_analyze(single([h_gt_4ch_cmr.PixelSpacing h_gt_4ch_cmr.SliceThickness]), CreateGtImageHeader([h_gt_4ch_cmr.PixelSpacing h_gt_4ch_cmr.SliceThickness]), fullfile(dst_dir, 'ch4_pixel_spacing'));
    Matlab_gt_write_analyze(single([1 1 h_gt_4ch_cmr.SliceThickness]), CreateGtImageHeader([1 1 h_gt_4ch_cmr.SliceThickness]), fullfile(dst_dir, 'ch4_1mm_pixel_spacing'));
    
    writeNPY(single(h_gt_2ch_cmr.ImagePositionPatient), fullfile(dst_dir, 'ch2_ImagePosition.npy'));
    writeNPY(single(h_gt_2ch_cmr.ImageOrientationPatient), fullfile(dst_dir, 'ch2_ImageOrientation.npy'));
    writeNPY(single([h_gt_2ch_cmr.PixelSpacing h_gt_2ch_cmr.SliceThickness]), fullfile(dst_dir, 'ch2_pixel_spacing.npy'));
    writeNPY(single([1 1 h_gt_2ch_cmr.SliceThickness]), fullfile(dst_dir, 'ch2_1mm_pixel_spacing.npy'));
    
    Matlab_gt_write_analyze(single(h_gt_2ch_cmr.ImagePositionPatient), CreateGtImageHeader(h_gt_2ch_cmr.ImagePositionPatient), fullfile(dst_dir, 'ch2_ImagePosition'));
    Matlab_gt_write_analyze(single(h_gt_2ch_cmr.ImageOrientationPatient), CreateGtImageHeader(h_gt_2ch_cmr.ImageOrientationPatient), fullfile(dst_dir, 'ch2_ImageOrientation'));
    Matlab_gt_write_analyze(single([h_gt_2ch_cmr.PixelSpacing h_gt_2ch_cmr.SliceThickness]), CreateGtImageHeader([h_gt_2ch_cmr.PixelSpacing h_gt_2ch_cmr.SliceThickness]), fullfile(dst_dir, 'ch2_pixel_spacing'));
    Matlab_gt_write_analyze(single([1 1 h_gt_2ch_cmr.SliceThickness]), CreateGtImageHeader([1 1 h_gt_2ch_cmr.SliceThickness]), fullfile(dst_dir, 'ch2_1mm_pixel_spacing'));
    
    sax_ImagePositionPatient = zeros(SLC, 3);
    for slc=1:SLC
        sax_ImagePositionPatient(slc,:) = h_gt_sax_cmr(slc).ImagePositionPatient;
    end
    
    writeNPY(single(sax_ImagePositionPatient), fullfile(dst_dir, 'sax_ImagePosition.npy'));
    writeNPY(single([h_gt_sax_cmr(1).PixelSpacing h_gt_sax_cmr(1).SliceThickness]), fullfile(dst_dir, 'sax_pixel_spacing.npy'));
    writeNPY(single(h_gt_sax_cmr(1).ImageOrientationPatient), fullfile(dst_dir, 'sax_ImageOrientation.npy'));
    writeNPY(single([1 1 h_gt_sax_cmr(1).SliceThickness]), fullfile(dst_dir, 'sax_1mm_pixel_spacing.npy'));
    
    Matlab_gt_write_analyze(single(sax_ImagePositionPatient), CreateGtImageHeader(sax_ImagePositionPatient), fullfile(dst_dir, 'sax_ImagePosition'));
    Matlab_gt_write_analyze(single([h_gt_sax_cmr(1).PixelSpacing h_gt_sax_cmr(1).SliceThickness]), CreateGtImageHeader([h_gt_sax_cmr(1).PixelSpacing h_gt_sax_cmr(1).SliceThickness]), fullfile(dst_dir, 'sax_pixel_spacing'));
    Matlab_gt_write_analyze(single(h_gt_sax_cmr(1).ImageOrientationPatient), CreateGtImageHeader(h_gt_sax_cmr(1).ImageOrientationPatient), fullfile(dst_dir, 'sax_ImageOrientation'));
    Matlab_gt_write_analyze(single([1 1 h_gt_sax_cmr(1).SliceThickness]), CreateGtImageHeader([1 1 h_gt_sax_cmr(1).SliceThickness]), fullfile(dst_dir, 'sax_1mm_pixel_spacing'));    
end
