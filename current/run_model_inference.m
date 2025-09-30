function model_list = run_model_inference(run_file_name, case_lists, gpu_lists)
% run_model_inference('/data/raw_data/SNR_PAPER/run_phantom.sh', case_lists)

model_list = {
    %"hrnet_more_T",             "/data/models/paper/mri_denoising_hrnet__whitney15_6xG8_MI300X_epoch_256_95p/whitney15_hrnet__more_T_6xG8_MI300X_256_95p_B4_Cine/checkpoints/mri-imaging-paper-STCNNT_HRNET_T1T1L1G1T1_T1T1L1G1T1_whitney15_hrnet__more_T_6xG8_MI300X_256_95p_B4_Cine_6x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_255.pth"; 
    "hrnet_more_T",             "/home/xueh/anaconda3/envs/gadgetron-gtprep/share/gadgetron/python/cmr_ml/models/hrnet_more_T__4xG16_V100_255_90p_B1_Cine_4x32G16-V100.pth"; 

    "hrnet",                    "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__baseline_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_HRNET_T1L1G1_T1L1G1_kings01_hrnet__baseline_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_59.pth"; 
    "hrnet_no_gmap",            "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__no_gmap_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_HRNET_T1L1G1_T1L1G1_kings01_hrnet__no_gmap_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_59.pth"; 
    "hrnet_no_MR_noise",        "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__no_MR_noise_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_HRNET_T1L1G1_T1L1G1_kings01_hrnet__no_MR_noise_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_59.pth"; 
    "hrnet_no_snr_unit",        "/data/models/paper/mri_denoising_hrnet__whitney15_2xG8_MI300X_epoch_60_95p/whitney15_hrnet__scale_by_signal_2xG8_MI300X_60_95p_B8_/checkpoints/mri-imaging-paper-STCNNT_HRNET_T1L1G1_T1L1G1_whitney15_hrnet__scale_by_signal_2xG8_MI300X_60_95p_B8__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_no_SNR_Unit/checkpoint_epoch_59.pth"; 
    "hrnet_dicom_mode",         "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__baseline_dicom_mode_2xG16_V100_60_95p_B1_/checkpoints/mri-imaging-paper-STCNNT_HRNET_T1L1G1_T1L1G1_kings01_hrnet__baseline_dicom_mode_2xG16_V100_60_95p_B1__2x32G16-V100_NN_32.0_32.0_C-64_residual_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_59.pth"; 
    
    "hrnet_TTT",                "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__TTT_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_HRNET_T1T1T1_T1T1T1_kings01_hrnet__TTT_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_59.pth"; 
    "hrnet_TTT_no_gmap",            "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__no_gmap_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_HRNET_T1L1G1_T1L1G1_kings01_hrnet__no_gmap_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_59.pth"; 
    "hrnet_TTT_no_MR_noise",        "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__no_MR_noise_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_HRNET_T1L1G1_T1L1G1_kings01_hrnet__no_MR_noise_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_59.pth"; 
    "hrnet_TTT_dicom_mode",         "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__baseline_dicom_mode_2xG16_V100_60_95p_B1_/checkpoints/mri-imaging-paper-STCNNT_HRNET_T1L1G1_T1L1G1_kings01_hrnet__baseline_dicom_mode_2xG16_V100_60_95p_B1__2x32G16-V100_NN_32.0_32.0_C-64_residual_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_59.pth"; 

    "hrnet_C3",                 "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__C3_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_HRNET_C3C3C3_C3C3C3_kings01_hrnet__C3_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_59.pth"; 
    "hrnet_C3_no_gmap",         "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__C3_no_gmap_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_HRNET_C3C3C3_C3C3C3_kings01_hrnet__C3_no_gmap_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_59.pth"; 
    "hrnet_C3_no_MR_noise",     "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__C3_no_MR_noise_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_HRNET_C3C3C3_C3C3C3_kings01_hrnet__C3_no_MR_noise_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_59.pth"; 
    "hrnet_C3_no_snr_unit",     "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__C3_scale_by_signal_2xG16_V100_60_95p_B4_/checkpoints/mri-imaging-paper-STCNNT_HRNET_C3C3C3_C3C3C3_kings01_hrnet__C3_scale_by_signal_2xG16_V100_60_95p_B4__2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_no_SNR_Unit/checkpoint_epoch_59.pth"; 
    "hrnet_C3_dicom_mode",      "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__C3_dicom_mode_2xG16_V100_60_95p_B1_/checkpoints/mri-imaging-paper-STCNNT_HRNET_C3C3C3_C3C3C3_kings01_hrnet__C3_dicom_mode_2xG16_V100_60_95p_B1__2x32G16-V100_NN_32.0_32.0_C-64_residual_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_59.pth"; 
    
    "hrnet_C2",                 "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__C2_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_HRNET_C2C2C2_C2C2C2_kings01_hrnet__C2_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_59.pth"; 
    
    "hrnet_ViT3D",              "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__ViT3D_2xG16_V100_60_95p_B1/checkpoints/mri-imaging-paper-STCNNT_HRNET_V3V3V3_V3V3V3_kings01_hrnet__ViT3D_2xG16_V100_60_95p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_59.pth"; 
    "hrnet_ViT3D_no_gmap",      "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__ViT3D_no_gmap_2xG16_V100_60_95p_B1/checkpoints/mri-imaging-paper-STCNNT_HRNET_V3V3V3_V3V3V3_kings01_hrnet__ViT3D_no_gmap_2xG16_V100_60_95p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_59.pth"; 
    "hrnet_ViT3D_no_MR_noise",  "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__ViT3D_no_MR_noise_2xG16_V100_60_95p_B1/checkpoints/mri-imaging-paper-STCNNT_HRNET_V3V3V3_V3V3V3_kings01_hrnet__ViT3D_no_MR_noise_2xG16_V100_60_95p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_59.pth"; 
    "hrnet_ViT3D_no_snr_unit",  "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__ViT3D_scale_by_signal_2xG16_V100_60_95p_B4_/checkpoints/mri-imaging-paper-STCNNT_HRNET_V3V3V3_V3V3V3_kings01_hrnet__ViT3D_scale_by_signal_2xG16_V100_60_95p_B4__2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_no_SNR_Unit/checkpoint_epoch_59.pth"; 
    "hrnet_ViT3D_dicom_mode",   "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__ViT3D_dicom_mode_2xG16_V100_60_95p_B1_/checkpoints/mri-imaging-paper-STCNNT_HRNET_V3V3V3_V3V3V3_kings01_hrnet__ViT3D_dicom_mode_2xG16_V100_60_95p_B1__2x32G16-V100_NN_32.0_32.0_C-64_residual_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_59.pth"; 
    
    "hrnet_ViT2D",              "/data/models/paper/mri_denoising_hrnet__kings01_2xG16_V100_epoch_60_95p/kings01_hrnet__ViT2D_2xG16_V100_60_95p_B1/checkpoints/mri-imaging-paper-STCNNT_HRNET_V2V2V2_V2V2V2_kings01_hrnet__ViT2D_2xG16_V100_60_95p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_59.pth"; 
    
    "unet",                     "/data/models/paper/mri_denoising_unet__kings01_2xG16_V100_epoch_60_95p/kings01_unet__baseline_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_UNET_T1L1G1_T1L1G1_kings01_unet__baseline_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_59.pth"; 
    "unet_no_gmap",             "/data/models/paper/mri_denoising_unet__kings01_2xG16_V100_epoch_60_95p/kings01_unet__no_gmap_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_UNET_T1L1G1_T1L1G1_kings01_unet__no_gmap_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_59.pth"; 
    "unet_no_MR_noise",         "/data/models/paper/mri_denoising_unet__kings01_2xG16_V100_epoch_60_95p/kings01_unet__no_MR_noise_2xG16_V100_60_95p_B1_Cine/checkpoints/mri-imaging-paper-STCNNT_UNET_T1L1G1_T1L1G1_kings01_unet__no_MR_noise_2xG16_V100_60_95p_B1_Cine_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_59.pth"; 
    "unet_no_snr_unit",         "/data/models/paper/mri_denoising_unet__whitney15_2xG8_MI300X_epoch_60_95p/whitney15_unet__scale_by_signal_2xG8_MI300X_60_95p_B8_/checkpoints/mri-imaging-paper-STCNNT_UNET_T1L1G1_T1L1G1_whitney15_unet__scale_by_signal_2xG8_MI300X_60_95p_B8__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_no_SNR_Unit/checkpoint_epoch_59.pth"; 
    "unet_dicom_mode",          "/data/models/paper/mri_denoising_unet__kings01_2xG16_V100_epoch_60_95p/kings01_unet__baseline_dicom_mode_2xG16_V100_60_95p_B1_/checkpoints/mri-imaging-paper-STCNNT_UNET_T1L1G1_T1L1G1_kings01_unet__baseline_dicom_mode_2xG16_V100_60_95p_B1__2x32G16-V100_NN_32.0_32.0_C-64_residual_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_59.pth"; 
    
    "unet_TTT",                 "/data/models/paper/mri_denoising_unet_whitney15_2xG8_MI300X_epoch_60_95p/whitney15_unet_TTT_2xG8_MI300X_60_95p_B2/checkpoints/mri-imaging-paper-STCNNT_UNET_T1T1T1_T1T1T1_whitney15_unet_TTT_2xG8_MI300X_60_95p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_59.pth";

    "unet_C3",                  "/data/models/paper/mri_denoising_unet_kings01_2xG16_V100_epoch_60_95p/kings01_unet_C3_2xG16_V100_60_95p_B1/checkpoints/mri-imaging-paper-STCNNT_UNET_C3C3C3_C3C3C3_kings01_unet_C3_2xG16_V100_60_95p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_59.pth"; 
    "unet_C3_no_gmap",          "/data/models/paper/mri_denoising_unet_whitney15_2xG8_MI300X_epoch_60_95p/whitney15_unet_C3_no_gmap_2xG8_MI300X_60_95p_B2/checkpoints/mri-imaging-paper-STCNNT_UNET_C3C3C3_C3C3C3_whitney15_unet_C3_no_gmap_2xG8_MI300X_60_95p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_59.pth"; 
    "unet_C3_no_MR_noise",      "/data/models/paper/mri_denoising_unet_whitney15_2xG8_MI300X_epoch_60_95p/whitney15_unet_C3_no_MR_noise_2xG8_MI300X_60_95p_B2/checkpoints/mri-imaging-paper-STCNNT_UNET_C3C3C3_C3C3C3_whitney15_unet_C3_no_MR_noise_2xG8_MI300X_60_95p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_59.pth"; 
    "unet_C3_no_snr_unit",      "/data/models/paper/mri_denoising_unet_kings01_2xG16_V100_epoch_60_95p/kings01_unet_C3_scale_by_signal_2xG16_V100_60_95p_B4_/checkpoints/mri-imaging-paper-STCNNT_UNET_C3C3C3_C3C3C3_kings01_unet_C3_scale_by_signal_2xG16_V100_60_95p_B4__2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_no_SNR_Unit/checkpoint_epoch_59.pth"; 
    "unet_C3_dicom_mode",       "/data/models/paper/mri_denoising_unet__kings01_2xG16_V100_epoch_60_95p/kings01_unet__C3_dicom_mode_2xG16_V100_60_95p_B1_/checkpoints/mri-imaging-paper-STCNNT_UNET_C3C3C3_C3C3C3_kings01_unet__C3_dicom_mode_2xG16_V100_60_95p_B1__2x32G16-V100_NN_32.0_32.0_C-64_residual_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_59.pth"; 
    
    "unet_C2",                  "/data/models/paper/mri_denoising_unet_whitney15_2xG8_MI300X_epoch_60_95p/whitney15_unet_C2_2xG8_MI300X_60_95p_B2/checkpoints/mri-imaging-paper-STCNNT_UNET_C2C2C2_C2C2C2_whitney15_unet_C2_2xG8_MI300X_60_95p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_59.pth"; 
    
    "unet_ViT3D",               "/data/models/paper/mri_denoising_unet__kings01_2xG16_V100_epoch_60_95p/kings01_unet__ViT3D_2xG16_V100_60_95p_B1/checkpoints/mri-imaging-paper-STCNNT_UNET_V3V3V3_V3V3V3_kings01_unet__ViT3D_2xG16_V100_60_95p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_59.pth"; 
    "unet_ViT3D_no_gmap",       "/data/models/paper/mri_denoising_unet__kings01_2xG16_V100_epoch_60_95p/kings01_unet__ViT3D_no_gmap_2xG16_V100_60_95p_B1/checkpoints/mri-imaging-paper-STCNNT_UNET_V3V3V3_V3V3V3_kings01_unet__ViT3D_no_gmap_2xG16_V100_60_95p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_59.pth"; 
    "unet_ViT3D_no_MR_noise",   "/data/models/paper/mri_denoising_unet__kings01_2xG16_V100_epoch_60_95p/kings01_unet__ViT3D_no_MR_noise_2xG16_V100_60_95p_B1/checkpoints/mri-imaging-paper-STCNNT_UNET_V3V3V3_V3V3V3_kings01_unet__ViT3D_no_MR_noise_2xG16_V100_60_95p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_59.pth"; 
    "unet_ViT3D_no_snr_unit",   "/data/models/paper/mri_denoising_unet__kings01_2xG16_V100_epoch_60_95p/kings01_unet__ViT3D_scale_by_signal_2xG16_V100_60_95p_B1/checkpoints/mri-imaging-paper-STCNNT_UNET_V3V3V3_V3V3V3_kings01_unet__ViT3D_scale_by_signal_2xG16_V100_60_95p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_no_SNR_Unit/checkpoint_epoch_59.pth"; 
    "unet_ViT3D_dicom_mode",    "/data/models/paper/mri_denoising_unet__kings01_2xG16_V100_epoch_60_95p/kings01_unet__ViT3D_dicom_mode_2xG16_V100_60_95p_B1__/checkpoints/mri-imaging-paper-STCNNT_UNET_V3V3V3_V3V3V3_kings01_unet__ViT3D_dicom_mode_2xG16_V100_60_95p_B1___2x32G16-V100_NN_32.0_32.0_C-64_residual_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_59.pth"; 
    
    "unet_ViT2D",               ""; 
    };

[fpath, fname, ext] = fileparts(run_file_name);

for g=1:numel(gpu_lists)
    gpu_str = num2str(gpu_lists(g));
    a_run_name = fullfile(fpath, [fname '_' gpu_str ext ]);

    f = fopen(a_run_name, "w")
    fprintf(f, '%s\n', "#!/bin/bash");
    fprintf(f, '%s\n', strjoin(["export CUDA_VISIBLE_DEVICES=" gpu_str], ""));
    fprintf(f, '%s\n', "export DISABLE_FLOAT16_INFERENCE=True");
    fprintf(f, '%s\n', "export RES_DIR=model_res");
    fprintf(f, '%s\n', 'echo "------------------------------------------------------------------------"');
    
    for i=g:numel(gpu_lists):numel(case_lists)
        a_case = case_lists{i}
        for k=1:size(model_list, 1)
            a_model_str = model_list{k, 1};
            a_model = model_list{k, 2};
            if (strlength(a_model)<5)
                continue
            end

            if isempty(strfind(a_model_str, 'no_gmap'))
                no_gmap = 0;
            else
                no_gmap = 1;
            end

            if isempty(strfind(a_model_str, 'no_MR_noise'))
                no_MR_noise = 0;
            else
                no_MR_noise = 1;
            end

            if isempty(strfind(a_model_str, 'no_snr_unit'))
                no_snr_unit = 0;
            else
                no_snr_unit = 1;
            end

            if isempty(strfind(a_model_str, 'dicom_mode'))
                dicom_mode = 0;
            else
                dicom_mode = 1;
            end

            cmd_str = create_run_command(a_case, a_model_str, a_model, no_gmap, no_MR_noise, no_snr_unit, dicom_mode);
            fprintf(f, [cmd_str{1} '\n\n']);
        end
        fprintf(f, '\n#############################################\n\n');
    end
    
    fclose(f);
    
    disp(['run file is - ' a_run_name])
end

end

function cmd_str = create_run_command(a_case, a_model_str, a_model, no_gmap, no_MR_noise, no_snr_unit, dicom_mode)
    cmd_str = ["python3 ./projects/mri_imaging/inference/run_inference.py --input_dir "  a_case " --output_dir " fullfile(a_case, a_model_str) " --scaling_factor 1.0 --im_scaling 1.0 --gmap_scaling 1.0 --input_fname input --gmap_fname gmap --saved_model_path " a_model " --batch_size 8 "];
    if (no_snr_unit)
        cmd_str = [cmd_str " --scale_by_signal"];
    end
    if (no_gmap)
        cmd_str = [cmd_str " --no_gmap"];
    end
    if (dicom_mode)
        cmd_str = [cmd_str " --no_gmap --scale_by_signal"];
    end
    if ~isempty(strfind(a_model_str, "ViT"))
        cmd_str = [cmd_str " --pad_time"];
    end
    cmd_str = strjoin(cmd_str, " ")
end