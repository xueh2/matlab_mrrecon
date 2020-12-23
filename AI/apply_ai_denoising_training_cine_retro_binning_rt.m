function apply_ai_denoising_training_cine_retro_binning_rt(dataDir, resDir, aiDir, files_all, script_name, suffix, denoise_py, noise_std, model_file, checkprocessed)
% apply_ai_denoising_training_cine_retro_binning_rt(dataDir, resDir, aiDir, files_all, script_name, suffix, denoise_py, noise_std, model_file, checkprocessed)

if(isunix())
    ext = '.sh';
else
    ext = '.bat';
end

fid = fopen(fullfile(aiDir, [script_name '_denoise' suffix ext]), 'a+');
set_env_apply_ai(fid);

try
    for n = 1:size(files_all,1)

        [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(files_all{n});

        case_dir = fullfile(resDir, study_dates, files_all{n});    
        dst_dir = fullfile(aiDir, study_dates, files_all{n});

        if(exist(case_dir))        

            disp(['--> process ' num2str(n) ' out of ' num2str(size(files_all,1)) ' - ' files_all{n}]);

            if(exist(fullfile(dst_dir, 'im_real.npy')))
                if(~exist(fullfile(dst_dir, 'im.npy')))
                    im_real = readNPY(fullfile(dst_dir, 'im_real.npy'));
                    im_imag = readNPY(fullfile(dst_dir, 'im_imag.npy'));
                    im = complex(im_real, im_imag);
                    im = abs(im);

                    writeNPY(single(im), fullfile(dst_dir, 'im.npy'));
                end
                
                res_case_Dir = fullfile(dst_dir, 'res_ai');
                if(~exist(res_case_Dir))
                    mkdir(res_case_Dir);
                end
                % lax_scripts
                if(isunix())
                    command = ['python3 '];
                else
                    command = ['python '];
                end

                if(exist(fullfile(dst_dir, 'im.npy')))
                    command = [command denoise_py ' --input ' fullfile(dst_dir, 'im.npy') ... 
                        ' --output ' fullfile(res_case_Dir, ['denoised' suffix '.npy']) ' --gmap ' fullfile(dst_dir, 'gfactor.npy') ...
                        ' --model ' model_file ' --ut_mode 0 --noise_std ' num2str(noise_std) ];

                    if(~checkprocessed || ~exist(fullfile(res_case_Dir, ['denoised' suffix '.npy'])))
                        command
                        fprintf(fid, '%s\n', command);
                    end
                end
            end        
        end    
    end
catch
    fclose(fid);
end
fclose(fid);
