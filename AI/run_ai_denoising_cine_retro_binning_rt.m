function run_ai_denoising_cine_retro_binning_rt(aiDir, resDir, files_all, batch_file, model_file, scale_factors, use_gmap, use_complex, is_CNNT)
% run_ai_denoising_cine_retro_binning_rt(aiDir, resDir, files_all, batch_file, model_file, scale_factors, use_gmap, use_complex, is_CNNT)

mkdir(resDir)

inference_call = '/home/xueh2/mrprogs/gadgetron_CMR_ML-source/projects/denoising/inference/run_inference_case.py';
inference_cnnt_call = '/home/xueh2/mrprogs/CNNT/run_cnnt_denoising_inference.py';

fid = fopen(batch_file, 'w');

for n = 1:size(files_all,1)
    
    fname = files_all{n, 1};
    if(iscell(fname))
        fname = fname{1};
    end
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(fname);
    
    ai_cast_dir = fullfile(aiDir, study_dates, fname);
    
    ai_dst_dir = fullfile(resDir, study_dates, fname);
    mkdir(ai_dst_dir)
    
    if(is_CNNT)
        command = ['python3 ' inference_cnnt_call ...
            ' --noisy_image ' ai_cast_dir '/im' ...
            ' --gmap ' ai_cast_dir '/gfactor' ...
            ' --output '  ai_dst_dir '/output.h5 ' ... 
            ' --model_path ' model_file ...
            ' --scale_factor ' num2str(scale_factors), ...
            ' --device cuda --iter_processing_times 2'];
    else
        command = ['python3 ' inference_call ...
            ' --case ' ai_cast_dir ...
            ' --output '  ai_dst_dir '/output.h5 ' ... 
            ' --model ' model_file ...
            ' --scale_factor ' num2str(scale_factors) ...
            ' --use_gmap ' num2str(use_gmap) ...
            ' --use_complex ' num2str(use_complex)];
    end
    
    fprintf(fid, '\n%s\n', command);
    
end

fclose(fid);
