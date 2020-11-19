function prepare_ai_training_fat_water_multi_slice(dataDir, resDir, aiDir, files_all, fat_series_num, water_series_num)
% prepare_ai_training_fat_water_multi_slice(dataDir, resDir, aiDir, files_all, fat_series_num, water_series_num)

mkdir(aiDir)

pic_dir = fullfile(aiDir, 'jpg_pics');
mkdir(pic_dir)

for n = 1:size(files_all,1)
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(files_all{n});
    
    case_dir = fullfile(resDir, study_dates, files_all{n});    
    dst_dir = fullfile(aiDir, study_dates, files_all{n});
    
    if(exist(case_dir))        
    
        disp(['--> process ' num2str(n) ' out of ' num2str(size(files_all,1)) ' - ' files_all{n}]);
        
        if(exist(fullfile(dst_dir, 'water.hdr')))
            continue;
        end
        
        try
            [water, w_header, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, water_series_num, 1);
            [fat, f_header, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, fat_series_num, 1);
        catch
            continue;                
        end
        water = squeeze(water);
        fat = squeeze(fat);

        mkdir(dst_dir);

        writeNPY(single(water), fullfile(dst_dir, 'water.npy'));
        writeNPY(single(fat), fullfile(dst_dir, 'fat.npy'));
        save(fullfile(dst_dir, 'headers.mat'), 'w_header', 'f_header');
        
        mkdir(fullfile(dst_dir, 'fat'));
        mkdir(fullfile(dst_dir, 'water'));    

        RO = size(fat, 1);
        E1 = size(fat, 2);

        fat = Matlab_gt_resize_2D_image(double(fat), 2*RO, 2*E1, 5);
        water = Matlab_gt_resize_2D_image(double(water), 2*RO, 2*E1, 5);

        header = CreateGtImageHeader(fat);
        Matlab_gt_write_analyze(single(fat), header, fullfile(dst_dir, 'fat'));
        Matlab_gt_write_analyze(single(water), header, fullfile(dst_dir, 'water'));

        N = size(fat, 3);

        fat2 = fat / (0.25 * max(fat(:)));
        for ii=1:N
            im = fat2(:,:,ii);
            im = im';
            im = flipdim(im, 1);
            imwrite(uint8(im*255), fullfile(dst_dir, 'fat', ['fat_' num2str(ii) '.jpg']));
        end

        water2 = water / (0.5 * max(water(:)));
        for ii=1:N
            im = water2(:,:,ii);
            im = im';
            im = flipdim(im, 1);
            imwrite(uint8(im*255), fullfile(dst_dir, 'water', ['water_' num2str(ii) '.jpg']));
        end
    end    
end