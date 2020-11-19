function prepare_ai_denoising_training_cine_retro_binning_rt(dataDir, resDir, aiDir, files_all, im_series_num, gfactor_series_num)
% prepare_ai_denoising_training_cine_retro_binning_rt(dataDir, resDir, aiDir, files_all, im_series_num, gfactor_series_num)

mkdir(aiDir)

pic_dir = fullfile(aiDir, 'jpg_pics');
mkdir(pic_dir)

for n = 1:size(files_all,1)
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(files_all{n});
    
    case_dir = fullfile(resDir, study_dates, files_all{n});    
    dst_dir = fullfile(aiDir, study_dates, files_all{n});
    
    if(exist(case_dir))        
    
        disp(['--> process ' num2str(n) ' out of ' num2str(size(files_all,1)) ' - ' files_all{n}]);
        
        if(exist(fullfile(pic_dir, [files_all{n} '.jpg'])))
            continue;
        end
        
        if(exist(fullfile(dst_dir, 'im.npy')))
            im = readNPY(fullfile(dst_dir, 'im.npy'));
            gfactor = readNPY(fullfile(dst_dir, 'gfactor.npy'));
        else
            try
                [im, im_header, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, im_series_num, 1);
                [gfactor, gmap_header, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, gfactor_series_num, 1);
            catch
                continue;                
            end
            im = squeeze(im);
            gfactor = squeeze(gfactor);

            mkdir(dst_dir);

            writeNPY(single(im), fullfile(dst_dir, 'im.npy'));
            writeNPY(single(gfactor), fullfile(dst_dir, 'gfactor.npy'));
            save(fullfile(dst_dir, 'headers.mat'), 'im_header', 'gmap_header');
        end
        
        im = uint8(255 * im / (0.5*max(im(:))) );
        
        if(numel(size(gfactor))==2)
            gfactor = repmat(gfactor, [1 1 size(im, 3)]);
        end
        
        if(numel(size(gfactor))==3 & numel(size(im))==4)
            gfactor = repmat(gfactor, [1 1 1 size(im, 4)]);
        end
        
        if(numel(size(im))==3)
            h = figure; imagescn(cat(4, im, gfactor), [], [1 2], [12], 3);
        else
            SLC = size(im, 3);
            h = figure; imagescn(cat(3, im, gfactor), [], [2 SLC], [24], 4);
        end
        saveas(h, fullfile(pic_dir, [files_all{n} '.jpg']), 'jpg');
        closeall
    end    
end