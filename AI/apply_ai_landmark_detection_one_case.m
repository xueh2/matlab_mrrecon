function [command_2D, command_3D, pts_2D, pts_3D] = apply_ai_landmark_detection_one_case(data_file, res_dir, batch_size, lax_model, three_D_model, RO, E1, RO_3D, E1_3D)
% [command_2D, command_3D, pts_2D, pts_3D] = apply_ai_landmark_detection_one_case(data_file, batch_size, lax_model, three_D_model, RO, E1, RO_3D, E1_3D)

    if(isunix())
        command = ['python3 '];
    else
        command = ['python '];
    end
    
    if(~exist(res_dir))
        mkdir(res_dir)
    end
    
    [path, sname, ext] = fileparts(data_file);
    
    pts_2D = fullfile(res_dir, ['AI_pts' '_' sname '.npy']);    
    command_2D = [command 'cmr_landmark_detection.py --input ' data_file ... 
        ' --output ' pts_2D ' --prob ' fullfile(res_dir, ['AI_probs' '_' sname '.npy']) ...
        ' --im_used ' fullfile(res_dir, ['AI_im_used' '_' sname '.npy']) ...
        ' --model ' lax_model ' --batch_size ' num2str(batch_size) ' --lax 1 --use_3D 0 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO) ' --E1 ' num2str(E1) ];

    pts_3D = fullfile(res_dir, ['AI_pts_3D' '_' sname '.npy']);
    command_3D = [command 'cmr_landmark_detection.py --input ' data_file ... 
        ' --output ' pts_3D ...
        ' --prob ' fullfile(res_dir, ['AI_probs_3D' '_' sname '.npy']) ...
        ' --im_used ' fullfile(res_dir, ['AI_im_used_3D' '_' sname '.npy']) ...
        ' --model_3D ' three_D_model ' --lax 1 --use_3D 1 --smooth_pts 0 --cli_mode 1 --fill_missed 1 --RO ' num2str(RO_3D) ' --E1 ' num2str(E1_3D) ];

