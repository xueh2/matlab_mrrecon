function create_run_script_resys(run_file_name, case_lists, res_dir, model_file, config_file, scaling_factor, batch_size, ignore_check_existed, cuda_str)
% create_run_script_resys(run_file_name, case_lists, res_dir, model_file, config_file, scaling_factor, batch_size, ignore_check_existed, cuda_str)
% python3 ./src/ifm/mri/denoising/run_inference.py --input_dir ${data_dir} --output_dir ${data_dir}/${res_dir} --saved_model_path $model_file --saved_config_path $config_file --input_fname input --gmap_fname gmap --scaling_factor ${scaling_factor} --batch_size 1

[fpath, fname, ext] = fileparts(run_file_name);

a_run_name = run_file_name;

f = fopen(a_run_name, "w")
fprintf(f, '%s\n', "#!/bin/bash");
fprintf(f, '%s', ["export CUDA_VISIBLE_DEVICES=" cuda_str])
fprintf(f, '\n%s\n', 'cd /home/xueh/mrprogs/resys-main/imaging-fm-projects')
fprintf(f, '\n%s\n', 'echo "------------------------------------------------------------------------"');

for i=1:size(case_lists, 1)
    a_case = case_lists{i, 1}
    if ~ignore_check_existed && exist(fullfile(a_case, 'res_dir', 'output.npy'))
        continue
    end
    cmd_str = create_run_command(a_case, fullfile(a_case, res_dir), model_file, config_file, scaling_factor, batch_size);
    fprintf(f, [cmd_str{1} '\n\n']);
    fprintf(f, '\n#############################################\n\n');
end

fclose(f);
    
disp(['run file is - ' a_run_name])

end

function cmd_str = create_run_command(a_case, res, model_file, config_file, scaling_factor, batch_size)
    cmd_str = ["python3 ./src/ifm/mri/denoising/run_inference.py --input_dir "  a_case " --output_dir " res " --saved_model_path " model_file " --saved_config_path " config_file " --scaling_factor " num2str(scaling_factor) " --batch_size " num2str(batch_size)];    
    cmd_str = strjoin(cmd_str, " ");
end