function prepare_ai_denoising_training_perfusion(dataDir, resDir, aiDir, files_all, im_series_num, gfactor_series_num)
% prepare_ai_denoising_training_perfusion(dataDir, resDir, aiDir, files_all, im_series_num, gfactor_series_num)

mkdir(aiDir)

pic_dir = fullfile(aiDir, 'jpg_pics');
mkdir(pic_dir)

for n = 1:size(files_all,1)
    
    fname = files_all{n, 1};
    fname = fname{1};
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(fname);
    
    case_dir = fullfile(resDir, study_dates, fname);    
    dst_dir = fullfile(aiDir, study_dates, fname);
    
    if(exist(case_dir))        
    
        disp(['--> process ' num2str(n) ' out of ' num2str(size(files_all,1)) ' - ' fname]);
        
        if(exist(fullfile(pic_dir, [fname '.jpg'])))
            continue;
        end
        
        if(exist(fullfile(dst_dir, 'im.npy')))
            if(exist(fullfile(dst_dir, 'im_real.npy')))
                im_real = readNPY(fullfile(dst_dir, 'im_real.npy'));
                im_imag = readNPY(fullfile(dst_dir, 'im_imag.npy'));
                im = complex(im_real, im_imag);
                im = abs(im);
            else
                im = readNPY(fullfile(dst_dir, 'im.npy'));
            end
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

            if(isreal(im))
                writeNPY(single(im), fullfile(dst_dir, 'im.npy'));
            else
                writeNPY(single(real(im)), fullfile(dst_dir, 'im_real.npy'));
                writeNPY(single(imag(im)), fullfile(dst_dir, 'im_imag.npy'));
            end
            writeNPY(single(gfactor), fullfile(dst_dir, 'gfactor.npy'));
            save(fullfile(dst_dir, 'headers.mat'), 'im_header', 'gmap_header');
            
            im = abs(im);
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