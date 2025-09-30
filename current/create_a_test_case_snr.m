function create_a_test_case_snr(run_file_name, case_lists, which_gmap, ignore_check_existed, output_dir)
% create_a_test_case_snr(run_file_name, case_lists)

[fpath, fname, ext] = fileparts(run_file_name);

a_run_name = run_file_name;

f = fopen(a_run_name, "w")
fprintf(f, '%s\n', "#!/bin/bash");
fprintf(f, '%s\n', 'echo "------------------------------------------------------------------------"');

for i=1:size(case_lists, 1)
    a_case = case_lists{i, 1}
    if ~ignore_check_existed && exist(fullfile(a_case, 'res', 'slc_0/gmap.npy'))
        continue
    end
    cmd_str = create_run_command(a_case, which_gmap, output_dir);
    fprintf(f, [cmd_str{1} '\n\n']);
    fprintf(f, '\n#############################################\n\n');
end

fclose(f);
    
disp(['run file is - ' a_run_name])

end

function cmd_str = create_run_command(a_case, which_gmap, output_dir)
    cmd_str = ["python3 /home/xueh/mrprogs/imagingfm_BTCHW/projects/mri_imaging/data/create_test_set_mri_snr_level.py --input_dir "  a_case " --output_dir " fullfile(a_case, output_dir) " --which_gmap " num2str(which_gmap) " --gmap_fname gmaps_aug"];    
    cmd_str = strjoin(cmd_str, " ");
end