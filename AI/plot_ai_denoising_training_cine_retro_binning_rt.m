function plot_ai_denoising_training_cine_retro_binning_rt(dataDir, resDir, aiDir, files_all, suffix)
% plot_ai_denoising_training_cine_retro_binning_rt(dataDir, resDir, aiDir, files_all, suffix)

pic_dir = fullfile(aiDir, 'jpgs');
mkdir(pic_dir);

for n = 1:size(files_all,1)

    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(files_all{n});

    case_dir = fullfile(resDir, study_dates, files_all{n});    
    dst_dir = fullfile(aiDir, study_dates, files_all{n});

    if(exist(case_dir))        

        disp(['--> process ' num2str(n) ' out of ' num2str(size(files_all,1)) ' - ' files_all{n}]);

        if(exist(fullfile(dst_dir, 'im_real.npy')))
            res_case_Dir = fullfile(dst_dir, 'res_ai');

            if(exist(fullfile(dst_dir, 'im.npy')) & exist(fullfile(res_case_Dir, ['denoised' suffix '.npy'])))
                im = readNPY(fullfile(dst_dir, 'im.npy'));
                output = readNPY(fullfile(res_case_Dir, ['denoised' suffix '.npy']));

                if(size(im,4)~=size(output,4))
                    im = permute(im, [1 2 4 3]);
                end

                SLC = size(im, 4);
                h = figure;
                if(SLC>1)
                    imagescn(cat(4, im, output, im-output), [0 mean(im(:))*3.5], [3, SLC], [], 3);
                else
                    imagescn(cat(4, im, output, im-output), [0 mean(im(:))*3.5], [1, 3], [], 3);
                end

                saveas(h, fullfile(pic_dir, [files_all{n} suffix '.jpg']), 'jpg');
                saveas(h, fullfile(pic_dir, [files_all{n} suffix '.fig']), 'fig');
            end
        end        
    end    
end
