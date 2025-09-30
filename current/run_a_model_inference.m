function run_a_model_inference(run_file_name, case_lists, gpu_lists, model_file, res_dir, check_processed)
% run_a_model_inference(run_file_name, case_lists, gpu_lists, model_file, check_processed)

model_dir = getenv('model_dir')
data_dir = getenv('data_dir')
output_dir = getenv('output_dir')

disp(['model_dir ' model_dir])
disp(['data_dir ' data_dir])
disp(['output_dir ' output_dir])

code_dir = getenv('code_dir');
if isempty(code_dir)
    code_dir = '/home/xueh/mrprogs/imagingfm_BTCHW';    
end
disp(['code_dir ' output_dir])

batch_size = getenv('batch_size')
disp(['batch_size ' batch_size])

if ~isempty(output_dir)
    mkdir(output_dir);
end

if nargin<5
    check_processed = 1;
end

[fpath, fname, ext] = fileparts(run_file_name);

for g=1:numel(gpu_lists)
    if isnumeric(gpu_lists(g))
        gpu_str = num2str(gpu_lists(g));
    else
        gpu_str = gpu_lists(g);
    end

    a_run_name = fullfile(fpath, [fname '_gpus_' num2str(g) ext ]);

    f = fopen(a_run_name, "w")
    fprintf(f, '%s\n', "#!/bin/bash");
    fprintf(f, '%s\n', strjoin(["export CUDA_VISIBLE_DEVICES=", gpu_str], ""));
    fprintf(f, '%s\n', "export DISABLE_FLOAT16_INFERENCE=True");
    fprintf(f, '%s\n', "export RES_DIR=model_res");
    fprintf(f, '%s\n', 'echo "------------------------------------------------------------------------"');
    
    for i=g:numel(gpu_lists):size(case_lists, 1)
        a_case = case_lists{i, 1}        
        if ~exist(a_case)
            a_case = fullfile(data_dir, a_case);
        end

        a_output_case_dir = fullfile(a_case, res_dir);

        if check_processed
            if exist(fullfile(a_output_case_dir, 'output_imag.npy')) & exist(fullfile(a_output_case_dir, 'output_real.npy'))
                continue
            end
        end

        cmd_str = create_run_command(code_dir, a_case, a_output_case_dir, model_file, batch_size, check_processed);
        fprintf(f, [cmd_str{1} '\n\n']);

        fprintf(f, '\n#############################################\n\n');
    end
    
    fclose(f);
    
    disp(['run file is - ' a_run_name])
end

end

function cmd_str = create_run_command(code_dir, a_case, a_output_case_dir, model_file, batch_size, check_processed)

    code_str = strjoin([code_dir "/projects/mri_imaging/inference/run_inference.py --check_processed " num2str(check_processed)], "");
    cmd_str = ["python3 " code_str " --input_dir "  a_case " --output_dir " a_output_case_dir " --scaling_factor 1.0 --im_scaling 1.0 --gmap_scaling 1.0 --input_fname input --gmap_fname gmap --saved_model_path " model_file " --batch_size " batch_size];

    cmd_str = strjoin(cmd_str, " ")
end