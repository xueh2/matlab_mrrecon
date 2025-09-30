function model_list = run_model_inference_v2(run_file_name, nnodes, case_lists, gpu_lists, check_processed, added_noise_sd, rep, input_model_list)
% model_list = run_model_inference_v2(run_file_name, nnodes, case_lists, gpu_lists, check_processed, added_noise_sd, rep, input_model_list)

model_dir = getenv('model_dir')
data_dir = getenv('data_dir')
output_dir = getenv('output_dir')
local_output_dir = getenv('local_output_dir')

disp(['model_dir ' model_dir])
disp(['data_dir ' data_dir])
disp(['output_dir ' output_dir])
disp(['local_output_dir ' local_output_dir])

if isempty(local_output_dir) | numel(local_output_dir)<2
    local_output_dir = output_dir;
end

code_dir = getenv('code_dir');
if isempty(code_dir)
    code_dir = '/home/xueh/mrprogs/imagingfm_BTCHW';    
end
disp(['code_dir ' code_dir])

batch_size = getenv('batch_size')
if isempty(batch_size)
    batch_size = 1;
end
disp(['batch_size ' batch_size])

if nargin<5
    check_processed = 1;
end

if nargin<6
    added_noise_sd = 0.1;
end

if nargin<7
    rep = 1;
end

if nargin<8
    input_model_list = [];
end

model_list = {

    "hrnet_TLTG",            "/data/models/paper_v3/model_search2_hrnet_2__whitney15_2xG8_MI300X_epoch_120_90p/whitney15_hrnet_2__deployed_T1L1T1G1_64_4_4_16_2xG8_MI300X_120_90p_B2_/checkpoints/paper-STCNNT_HRNET_T1L1T1G1_T1L1T1G1_whitney15_hrnet_2__deployed_T1L1T1G1_64_4_4_16_2xG8_MI300X_120_90p_B2__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth",

    "hrnet_TLG",             "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_120_90p/kings01_hrnet_2__deployed_T1L1G1_32_4_4_16_2xG16_V100_120_90p_B1_/checkpoints/paper-STCNNT_HRNET_T1L1G1_T1L1G1_kings01_hrnet_2__deployed_T1L1G1_32_4_4_16_2xG16_V100_120_90p_B1__2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth"; 
    "hrnet_TLG_no_gmap",     "/data/models/paper_v3/model_search2_hrnet_2__whitney15_1xG8_MI300X_epoch_80_90p/whitney15_hrnet_2__T1L1G1_64_4_4_16_no_gmap_1xG8_MI300X_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1L1G1_T1L1G1_whitney15_hrnet_2__T1L1G1_64_4_4_16_no_gmap_1xG8_MI300X_80_90p_B2_1x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "hrnet_TLG_no_MR_noise", "/data/models/paper_v3/model_search2_hrnet_2__whitney15_1xG8_MI300X_epoch_80_90p/whitney15_hrnet_2__T1L1G1_64_4_4_16_no_MR_noise_1xG8_MI300X_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1L1G1_T1L1G1_whitney15_hrnet_2__T1L1G1_64_4_4_16_no_MR_noise_1xG8_MI300X_80_90p_B2_1x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "hrnet_TLG_dicom_mode",  "/data/models/paper_v3/model_search2_hrnet_2__whitney15_1xG8_MI300X_epoch_80_90p/whitney15_hrnet_2__T1L1G1_64_4_4_16_dicom_mode_1xG8_MI300X_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1L1G1_T1L1G1_whitney15_hrnet_2__T1L1G1_64_4_4_16_dicom_mode_1xG8_MI300X_80_90p_B2_1x192G8-MI300X_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "hrnet_TLGTLG",             "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_120_90p/kings01_hrnet_2__deployed_T1L1T1G1T1_32_4_4_16_2xG16_V100_120_90p_B1_/checkpoints/paper-STCNNT_HRNET_T1L1T1G1T1_T1L1T1G1T1_kings01_hrnet_2__deployed_T1L1T1G1T1_32_4_4_16_2xG16_V100_120_90p_B1__2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth"; 
    "hrnet_TLGTLG_no_gmap",     "/data/models/paper_v3/model_search2_hrnet_2__whitney15_1xG8_MI300X_epoch_80_90p/whitney15_hrnet_2__T1L1G1T1L1G1_64_4_4_16_no_gmap_1xG8_MI300X_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1L1G1T1L1G1_T1L1G1T1L1G1_whitney15_hrnet_2__T1L1G1T1L1G1_64_4_4_16_no_gmap_1xG8_MI300X_80_90p_B2_1x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "hrnet_TLGTLG_no_MR_noise", "/data/models/paper_v3/model_search2_hrnet_2__whitney15_1xG8_MI300X_epoch_80_90p/whitney15_hrnet_2__T1L1G1T1L1G1_64_4_4_16_no_MR_noise_1xG8_MI300X_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1L1G1T1L1G1_T1L1G1T1L1G1_whitney15_hrnet_2__T1L1G1T1L1G1_64_4_4_16_no_MR_noise_1xG8_MI300X_80_90p_B2_1x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "hrnet_TLGTLG_dicom_mode",  ""; 

    "hrnet_TTT",                "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1_T1T1T1_kings01_hrnet_2__T1T1T1_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "hrnet_TTT_no_gmap",        "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1_64_no_gmap_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1_T1T1T1_kings01_hrnet_2__T1T1T1_64_no_gmap_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "hrnet_TTT_no_MR_noise",    "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1_64_no_MR_noise_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1_T1T1T1_kings01_hrnet_2__T1T1T1_64_no_MR_noise_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "hrnet_TTT_dicom_mode",     "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1_64_dicom_mode_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1_T1T1T1_kings01_hrnet_2__T1T1T1_64_dicom_mode_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "hrnet_TTTTTT",             "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_hrnet_2__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "hrnet_TTTTTT_no_gmap",     "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1T1T1T1_64_no_gmap_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_hrnet_2__T1T1T1T1T1T1_64_no_gmap_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "hrnet_TTTTTT_no_MR_noise", "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1T1T1T1_64_no_MR_noise_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_hrnet_2__T1T1T1T1T1T1_64_no_MR_noise_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "hrnet_TTTTTT_dicom_mode",  "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1T1T1T1_64_dicom_mode_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_hrnet_2__T1T1T1T1T1T1_64_dicom_mode_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "hrnet_Swin",              "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_2xG16_V100_80_90p_B1/checkpoints/paper-STCNNT_HRNET_S3ShS3ShS3Sh_S3ShS3ShS3Sh_kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_2xG16_V100_80_90p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "hrnet_Swin_no_gmap",      "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_no_gmap_2xG16_V100_80_90p_B1/checkpoints/paper-STCNNT_HRNET_S3ShS3ShS3Sh_S3ShS3ShS3Sh_kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_no_gmap_2xG16_V100_80_90p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "hrnet_Swin_no_MR_noise",  "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_no_MR_noise_2xG16_V100_80_90p_B1/checkpoints/paper-STCNNT_HRNET_S3ShS3ShS3Sh_S3ShS3ShS3Sh_kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_no_MR_noise_2xG16_V100_80_90p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "hrnet_Swin_dicom_mode",   "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_dicom_mode_2xG16_V100_80_90p_B1/checkpoints/paper-STCNNT_HRNET_S3ShS3ShS3Sh_S3ShS3ShS3Sh_kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_dicom_mode_2xG16_V100_80_90p_B1_2x32G16-V100_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "hrnet_C3",                 "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__C3C3C3_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_C3C3C3_C3C3C3_kings01_hrnet_2__C3C3C3_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "hrnet_C3_no_gmap",         "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__C3C3C3_64_no_gmap_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_C3C3C3_C3C3C3_kings01_hrnet_2__C3C3C3_64_no_gmap_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "hrnet_C3_no_MR_noise",     "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__C3C3C3_64_no_MR_noise_2xG16_V100_80_90p_B4_/checkpoints/paper-STCNNT_HRNET_C3C3C3_C3C3C3_kings01_hrnet_2__C3C3C3_64_no_MR_noise_2xG16_V100_80_90p_B4__2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "hrnet_C3_dicom_mode",      "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__C3C3C3_64_dicom_mode_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_C3C3C3_C3C3C3_kings01_hrnet_2__C3C3C3_64_dicom_mode_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "hrnet_C2",                 "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__C2C2C2_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_C2C2C2_C2C2C2_kings01_hrnet_2__C2C2C2_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "hrnet_ViT3D",              "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__V3V3V3_64_8_8_16_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_V3V3V3_V3V3V3_kings01_hrnet_2__V3V3V3_64_8_8_16_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "hrnet_ViT3D_no_gmap",      "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__V3V3V3_64_8_8_16_no_gmap_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_V3V3V3_V3V3V3_kings01_hrnet_2__V3V3V3_64_8_8_16_no_gmap_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "hrnet_ViT3D_no_MR_noise",  "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__V3V3V3_64_8_8_16_no_MR_noise_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_V3V3V3_V3V3V3_kings01_hrnet_2__V3V3V3_64_8_8_16_no_MR_noise_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "hrnet_ViT3D_dicom_mode",   "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__V3V3V3_64_8_8_16_dicom_mode_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_V3V3V3_V3V3V3_kings01_hrnet_2__V3V3V3_64_8_8_16_dicom_mode_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "hrnet_ViT2D",              "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__V2V2V2_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_V2V2V2_V2V2V2_kings01_hrnet_2__V2V2V2_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    % ----------------------------------------------------------------

    "unet_TLG",             "/data/models/paper_v3/model_search_unet__msroctobasicvc_2xG8_MI300X_epoch_80_90p/msroctobasicvc_unet__T1L1G1_64_4_4_16_T1L1G1_2xG8_MI300X_80_90p_B2/checkpoints/mri-imaging-para-search-STCNNT_UNET_T1L1G1_T1L1G1_msroctobasicvc_unet__T1L1G1_64_4_4_16_T1L1G1_2xG8_MI300X_80_90p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "unet_TLG_no_gmap",     "/data/models/paper_v3/model_search_unet__msroctobasicvc_2xG8_MI300X_epoch_80_90p/msroctobasicvc_unet__T1L1G1_64_4_4_16_no_gmap_T1L1G1_2xG8_MI300X_80_90p_B2/checkpoints/mri-imaging-para-search-STCNNT_UNET_T1L1G1_T1L1G1_msroctobasicvc_unet__T1L1G1_64_4_4_16_no_gmap_T1L1G1_2xG8_MI300X_80_90p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "unet_TLG_no_MR_noise", "/data/models/paper_v3/model_search_unet__msroctobasicvc_2xG8_MI300X_epoch_80_90p/msroctobasicvc_unet__T1L1G1_64_4_4_16_no_MR_noise_T1L1G1_2xG8_MI300X_80_90p_B2/checkpoints/mri-imaging-para-search-STCNNT_UNET_T1L1G1_T1L1G1_msroctobasicvc_unet__T1L1G1_64_4_4_16_no_MR_noise_T1L1G1_2xG8_MI300X_80_90p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "unet_TLG_dicom_mode",  "/data/models/paper_v3/model_search_unet__msroctobasicvc_2xG8_MI300X_epoch_80_90p/msroctobasicvc_unet__T1L1G1_64_4_4_16_dicom_mode_T1L1G1_2xG8_MI300X_80_90p_B2/checkpoints/mri-imaging-para-search-STCNNT_UNET_T1L1G1_T1L1G1_msroctobasicvc_unet__T1L1G1_64_4_4_16_dicom_mode_T1L1G1_2xG8_MI300X_80_90p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "unet_TLGTLG",             "/data/models/paper_v3/model_search_unet__msroctobasicvc_2xG8_MI300X_epoch_80_90p/msroctobasicvc_unet__T1L1G1T1L1G1_64_4_4_16_T1L1G1T1L1G1_2xG8_MI300X_80_90p_B2/checkpoints/mri-imaging-para-search-STCNNT_UNET_T1L1G1T1L1G1_T1L1G1T1L1G1_msroctobasicvc_unet__T1L1G1T1L1G1_64_4_4_16_T1L1G1T1L1G1_2xG8_MI300X_80_90p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "unet_TLGTLG_no_gmap",     "/data/models/paper_v3/model_search_unet__msroctobasicvc_2xG8_MI300X_epoch_80_90p/msroctobasicvc_unet__T1L1G1T1L1G1_64_4_4_16_no_gmap_T1L1G1T1L1G1_2xG8_MI300X_80_90p_B2/checkpoints/mri-imaging-para-search-STCNNT_UNET_T1L1G1T1L1G1_T1L1G1T1L1G1_msroctobasicvc_unet__T1L1G1T1L1G1_64_4_4_16_no_gmap_T1L1G1T1L1G1_2xG8_MI300X_80_90p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "unet_TLGTLG_no_MR_noise", "/data/models/paper_v3/model_search_unet__msroctobasicvc_2xG8_MI300X_epoch_80_90p/msroctobasicvc_unet__T1L1G1T1L1G1_64_4_4_16_no_MR_noise_T1L1G1T1L1G1_2xG8_MI300X_80_90p_B2/checkpoints/mri-imaging-para-search-STCNNT_UNET_T1L1G1T1L1G1_T1L1G1T1L1G1_msroctobasicvc_unet__T1L1G1T1L1G1_64_4_4_16_no_MR_noise_T1L1G1T1L1G1_2xG8_MI300X_80_90p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "unet_TLGTLG_dicom_mode",  "/data/models/paper_v3/model_search_unet__msroctobasicvc_2xG8_MI300X_epoch_80_90p/msroctobasicvc_unet__T1L1G1T1L1G1_64_4_4_16_dicom_mode_T1L1G1T1L1G1_2xG8_MI300X_80_90p_B2/checkpoints/mri-imaging-para-search-STCNNT_UNET_T1L1G1T1L1G1_T1L1G1T1L1G1_msroctobasicvc_unet__T1L1G1T1L1G1_64_4_4_16_dicom_mode_T1L1G1T1L1G1_2xG8_MI300X_80_90p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "unet_TTT",                "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__T1T1T1_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_T1T1T1_T1T1T1_kings01_unet__T1T1T1_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "unet_TTT_no_gmap",        "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__T1T1T1_64_no_gmap_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_T1T1T1_T1T1T1_kings01_unet__T1T1T1_64_no_gmap_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "unet_TTT_no_MR_noise",    "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__T1T1T1_64_no_MR_noise_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_T1T1T1_T1T1T1_kings01_unet__T1T1T1_64_no_MR_noise_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "unet_TTT_dicom_mode",     "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__T1T1T1_64_dicom_mode_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_T1T1T1_T1T1T1_kings01_unet__T1T1T1_64_dicom_mode_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "unet_TTTTTT",             "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_unet__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "unet_TTTTTT_no_gmap",     "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__T1T1T1T1T1T1_64_no_gmap_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_unet__T1T1T1T1T1T1_64_no_gmap_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "unet_TTTTTT_no_MR_noise", "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__T1T1T1T1T1T1_64_no_MR_noise_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_unet__T1T1T1T1T1T1_64_no_MR_noise_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "unet_TTTTTT_dicom_mode",  "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__T1T1T1T1T1T1_64_dicom_mode_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_unet__T1T1T1T1T1T1_64_dicom_mode_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "unet_Swin",              "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__S3ShS3ShS3Sh_64_8_8_16_2xG16_V100_80_90p_B1/checkpoints/paper-STCNNT_UNET_S3ShS3ShS3Sh_S3ShS3ShS3Sh_kings01_unet__S3ShS3ShS3Sh_64_8_8_16_2xG16_V100_80_90p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "unet_Swin_no_gmap",      "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__S3ShS3ShS3Sh_64_8_8_16_no_gmap_2xG16_V100_80_90p_B1/checkpoints/paper-STCNNT_UNET_S3ShS3ShS3Sh_S3ShS3ShS3Sh_kings01_unet__S3ShS3ShS3Sh_64_8_8_16_no_gmap_2xG16_V100_80_90p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "unet_Swin_no_MR_noise",  "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__S3ShS3ShS3Sh_64_8_8_16_no_MR_noise_2xG16_V100_80_90p_B1/checkpoints/paper-STCNNT_UNET_S3ShS3ShS3Sh_S3ShS3ShS3Sh_kings01_unet__S3ShS3ShS3Sh_64_8_8_16_no_MR_noise_2xG16_V100_80_90p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "unet_Swin_dicom_mode",   "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__S3ShS3ShS3Sh_64_8_8_16_dicom_mode_2xG16_V100_80_90p_B1/checkpoints/paper-STCNNT_UNET_S3ShS3ShS3Sh_S3ShS3ShS3Sh_kings01_unet__S3ShS3ShS3Sh_64_8_8_16_dicom_mode_2xG16_V100_80_90p_B1_2x32G16-V100_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "unet_C3",                 "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__C3C3C3_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_C3C3C3_C3C3C3_kings01_unet__C3C3C3_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "unet_C3_no_gmap",         "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__C3C3C3_64_no_gmap_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_C3C3C3_C3C3C3_kings01_unet__C3C3C3_64_no_gmap_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "unet_C3_no_MR_noise",     "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__C3C3C3_64_no_MR_noise_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_C3C3C3_C3C3C3_kings01_unet__C3C3C3_64_no_MR_noise_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "unet_C3_dicom_mode",      "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__C3C3C3_64_dicom_mode_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_C3C3C3_C3C3C3_kings01_unet__C3C3C3_64_dicom_mode_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 

    "unet_C2",                 "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__C2C2C2_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_C2C2C2_C2C2C2_kings01_unet__C2C2C2_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "unet_ViT3D",              "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__V3V3V3_64_8_8_16_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_V3V3V3_V3V3V3_kings01_unet__V3V3V3_64_8_8_16_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "unet_ViT3D_no_gmap",      "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__V3V3V3_64_8_8_16_no_gmap_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_V3V3V3_V3V3V3_kings01_unet__V3V3V3_64_8_8_16_no_gmap_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "unet_ViT3D_no_MR_noise",  "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__V3V3V3_64_8_8_16_no_MR_noise_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_V3V3V3_V3V3V3_kings01_unet__V3V3V3_64_8_8_16_no_MR_noise_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "unet_ViT3D_dicom_mode",   "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__V3V3V3_64_8_8_16_dicom_mode_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_V3V3V3_V3V3V3_kings01_unet__V3V3V3_64_8_8_16_dicom_mode_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 
    "unet_ViT3D_dicom_mode_2",    "/data/models/paper/mri_denoising_unet__kings01_2xG16_V100_epoch_60_95p/kings01_unet__ViT3D_dicom_mode_2xG16_V100_60_95p_B1__/checkpoints/mri-imaging-paper-STCNNT_UNET_V3V3V3_V3V3V3_kings01_unet__ViT3D_dicom_mode_2xG16_V100_60_95p_B1___2x32G16-V100_NN_32.0_32.0_C-64_residual_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_59.pth"; 

    "unet_ViT2D",              "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__V2V2V2_64_2xG16_V100_80_90p_B2_/checkpoints/paper-STCNNT_UNET_V2V2V2_V2V2V2_kings01_unet__V2V2V2_64_2xG16_V100_80_90p_B2__2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    };

model_list = {

    %"hrnet_TLGTLG",             "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_120_90p/kings01_hrnet_2__deployed_T1L1T1G1T1_32_4_4_16_2xG16_V100_120_90p_B1_/checkpoints/paper-STCNNT_HRNET_T1L1T1G1T1_T1L1T1G1T1_kings01_hrnet_2__deployed_T1L1T1G1T1_32_4_4_16_2xG16_V100_120_90p_B1__2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth"; 
    "hrnet_TLGTLG",             "/data/models/paper_v3/model_search2_hrnet_2__whitney15_2xG8_MI300X_epoch_120_90p/whitney15_hrnet_2__deployed_T1L1G1T1L1G1_64_4_4_8_2xG8_MI300X_120_90p_B2_/checkpoints/paper-STCNNT_HRNET_T1L1G1T1L1G1_T1L1G1T1L1G1_whitney15_hrnet_2__deployed_T1L1G1T1L1G1_64_4_4_8_2xG8_MI300X_120_90p_B2__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth"

    "hrnet_TTTTTT",             "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_hrnet_2__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "hrnet_TTT",                "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1_T1T1T1_kings01_hrnet_2__T1T1T1_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "hrnet_TLTG",            "/data/models/paper_v3/model_search2_hrnet_2__whitney15_2xG8_MI300X_epoch_120_90p/whitney15_hrnet_2__deployed_T1L1T1G1_64_4_4_16_2xG8_MI300X_120_90p_B2_/checkpoints/paper-STCNNT_HRNET_T1L1T1G1_T1L1T1G1_whitney15_hrnet_2__deployed_T1L1T1G1_64_4_4_16_2xG8_MI300X_120_90p_B2__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth",

    "hrnet_TLG",             "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_120_90p/kings01_hrnet_2__deployed_T1L1G1_32_4_4_16_2xG16_V100_120_90p_B1_/checkpoints/paper-STCNNT_HRNET_T1L1G1_T1L1G1_kings01_hrnet_2__deployed_T1L1G1_32_4_4_16_2xG16_V100_120_90p_B1__2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth"; 


    "hrnet_Swin",              "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_2xG16_V100_80_90p_B1/checkpoints/paper-STCNNT_HRNET_S3ShS3ShS3Sh_S3ShS3ShS3Sh_kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_2xG16_V100_80_90p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "hrnet_C3",                 "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__C3C3C3_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_C3C3C3_C3C3C3_kings01_hrnet_2__C3C3C3_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "hrnet_C2",                 "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__C2C2C2_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_C2C2C2_C2C2C2_kings01_hrnet_2__C2C2C2_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "hrnet_ViT3D",              "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__V3V3V3_64_8_8_16_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_V3V3V3_V3V3V3_kings01_hrnet_2__V3V3V3_64_8_8_16_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "hrnet_ViT2D",              "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__V2V2V2_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_V2V2V2_V2V2V2_kings01_hrnet_2__V2V2V2_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    %----------------------------------------------------------------

    %"unet_TLG",             "/data/models/paper_v3/model_search2_unet__whitney15_2xG8_MI300X_epoch_80_90p/whitney15_unet__T1L1G1_64_4_4_16_2xG8_MI300X_80_90p_B4/checkpoints/paper-STCNNT_UNET_T1L1G1_T1L1G1_whitney15_unet__T1L1G1_64_4_4_16_2xG8_MI300X_80_90p_B4_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual//checkpoint_epoch_79.pth"; 
    "unet_TLG",             "/data/models/paper_v3/model_search2_unet__whitney15_2xG8_MI300X_epoch_120_95p/whitney15_unet__deployed_T1L1G1_64_4_4_16_2xG8_MI300X_120_95p_B2_/checkpoints/paper-STCNNT_UNET_T1L1G1_T1L1G1_whitney15_unet__deployed_T1L1G1_64_4_4_16_2xG8_MI300X_120_95p_B2__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth";

    %"unet_TLGTLG",             "/data/models/paper_v3/model_search_unet__msroctobasicvc_2xG8_MI300X_epoch_80_90p/msroctobasicvc_unet__T1L1G1T1L1G1_64_4_4_16_T1L1G1T1L1G1_2xG8_MI300X_80_90p_B2/checkpoints/mri-imaging-para-search-STCNNT_UNET_T1L1G1T1L1G1_T1L1G1T1L1G1_msroctobasicvc_unet__T1L1G1T1L1G1_64_4_4_16_T1L1G1T1L1G1_2xG8_MI300X_80_90p_B2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual//checkpoint_epoch_79.pth"; 

    "unet_TLGTLG",             "/data/models/paper_v3/model_search2_unet__whitney15_2xG8_MI300X_epoch_128_95p/whitney15_unet__deployed_T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_128_95p_B2_/checkpoints/paper-STCNNT_UNET_T1L1G1T1L1G1_T1L1G1T1L1G1_whitney15_unet__deployed_T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_128_95p_B2__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_127.pth";

    "unet_TTT",                "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__T1T1T1_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_T1T1T1_T1T1T1_kings01_unet__T1T1T1_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "unet_TTTTTT",             "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_unet__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "unet_Swin",              "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__S3ShS3ShS3Sh_64_8_8_16_2xG16_V100_80_90p_B1/checkpoints/paper-STCNNT_UNET_S3ShS3ShS3Sh_S3ShS3ShS3Sh_kings01_unet__S3ShS3ShS3Sh_64_8_8_16_2xG16_V100_80_90p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "unet_C3",                 "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__C3C3C3_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_C3C3C3_C3C3C3_kings01_unet__C3C3C3_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "unet_C2",                 "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__C2C2C2_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_C2C2C2_C2C2C2_kings01_unet__C2C2C2_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "unet_ViT3D",              "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__V3V3V3_64_8_8_16_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_UNET_V3V3V3_V3V3V3_kings01_unet__V3V3V3_64_8_8_16_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    "unet_ViT2D",              "/data/models/paper_v3/model_search2_unet__kings01_2xG16_V100_epoch_80_90p/kings01_unet__V2V2V2_64_2xG16_V100_80_90p_B2_/checkpoints/paper-STCNNT_UNET_V2V2V2_V2V2V2_kings01_unet__V2V2V2_64_2xG16_V100_80_90p_B2__2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 

    };

model_list = {

    "hrnet_TLGTLG",             "/data/models/paper_v3/model_search2_hrnet_2__whitney15_2xG8_MI300X_epoch_120_90p/whitney15_hrnet_2__deployed_T1L1G1T1L1G1_64_4_4_8_2xG8_MI300X_120_90p_B2_/checkpoints/paper-STCNNT_HRNET_T1L1G1T1L1G1_T1L1G1T1L1G1_whitney15_hrnet_2__deployed_T1L1G1T1L1G1_64_4_4_8_2xG8_MI300X_120_90p_B2__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth"

    "hrnet_TLGTLGTLGTLG",       "/data/models/paper_v3/model_search2_hrnet_2__whitney15_2xG8_MI300X_epoch_120_95p/whitney15_hrnet_2__deployed_T1L1G1T1L1G1T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_120_95p_B2_/checkpoints/paper-STCNNT_HRNET_T1L1G1T1L1G1T1L1G1T1L1G1_T1L1G1T1L1G1T1L1G1T1L1G1_whitney15_hrnet_2__deployed_T1L1G1T1L1G1T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_120_95p_B2__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth"

    "hrnet_TLGTLGTLGTLG_SP0",   "/data/models/paper_v3/model_search2_hrnet__whitney15_2xG8_MI300X_epoch_160_95p/whitney15_hrnet__deployed_T1L1G1T1L1G1T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_160_95p_B2_layer_sp0/checkpoints/paper-STCNNT_HRNET_T1L1G1T1L1G1T1L1G1T1L1G1_T1L1G1T1L1G1T1L1G1T1L1G1_whitney15_hrnet__deployed_T1L1G1T1L1G1T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_160_95p_B2_layer_sp0_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/best_checkpoint_epoch_140.pth"

    "hrnet_TLGTLG",       "/data/models/paper_v3/model_search2_hrnet__whitney15_2xG8_MI300X_epoch_120_95p/whitney15_hrnet__deployed_with_2d_T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_120_95p_B2_layer_sp0/checkpoints/paper-STCNNT_HRNET_T1L1G1T1L1G1_T1L1G1T1L1G1_whitney15_hrnet__deployed_with_2d_T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_120_95p_B2_layer_sp0_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth"

    %"unet_TLGTLG",             "/data/models/paper_v3/model_search2_unet__whitney15_2xG8_MI300X_epoch_128_95p/whitney15_unet__deployed_T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_128_95p_B2_/checkpoints/paper-STCNNT_UNET_T1L1G1T1L1G1_T1L1G1T1L1G1_whitney15_unet__deployed_T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_128_95p_B2__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_127.pth";

    %----------------------------------------------------------------

    %"unet_TLG",             "/data/models/paper_v3/model_search2_unet__whitney15_2xG8_MI300X_epoch_120_95p/whitney15_unet__deployed_T1L1G1_64_4_4_16_2xG8_MI300X_120_95p_B2_/checkpoints/paper-STCNNT_UNET_T1L1G1_T1L1G1_whitney15_unet__deployed_T1L1G1_64_4_4_16_2xG8_MI300X_120_95p_B2__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth";

    %"unet_TLGTLG",             "/data/models/paper_v3/model_search2_hrnet_2__whitney15_2xG8_MI300X_epoch_128_98p/whitney15_hrnet_2__deployed_with_2d_T1L1T1G1_64_4_4_16_2xG8_MI300X_128_98p_B4_/checkpoints/paper-STCNNT_HRNET_T1L1T1G1_T1L1T1G1_whitney15_hrnet_2__deployed_with_2d_T1L1T1G1_64_4_4_16_2xG8_MI300X_128_98p_B4__2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_119.pth";
    };


model_list = {
    "hrnet_TTTTTT",             "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_hrnet_2__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "hrnet_TTTTTT_no_gmap",     "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1T1T1T1_64_no_gmap_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_hrnet_2__T1T1T1T1T1T1_64_no_gmap_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_ignore_gmap/checkpoint_epoch_79.pth"; 
    "hrnet_TTTTTT_no_MR_noise", "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1T1T1T1_64_no_MR_noise_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_hrnet_2__T1T1T1T1T1T1_64_no_MR_noise_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual_only_white_noise/checkpoint_epoch_79.pth"; 
    "hrnet_TTTTTT_dicom_mode",  "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1T1T1T1_64_dicom_mode_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_hrnet_2__T1T1T1T1T1T1_64_dicom_mode_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_residual_only_white_noise_ignore_gmap_no_SNR_Unit_dicom_mode/checkpoint_epoch_79.pth"; 
    };

model_list = {
    "hrnet_28m", "/data/models/paper_v3/model_search2_hrnet__whitney15_2xG8_MI300X_epoch_160_95p/whitney15_hrnet__deployed_T1L1G1_64_ws4_4_16_ps2_2_2_2xG8_MI300X_160_95p_B1_layer_sp0_patch_size2/checkpoints/paper-STCNNT_HRNET_T1L1G1_T1L1G1_whitney15_hrnet__deployed_T1L1G1_64_ws4_4_16_ps2_2_2_2xG8_MI300X_160_95p_B1_layer_sp0_patch_size2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_159.pth";
    "hrnet_50m",  "/data/models/paper_v3/model_search2_hrnet__whitney15_2xG8_MI300X_epoch_160_95p/whitney15_hrnet__deployed_T1L1G1T1L1G1_64_ws4_4_16_ps2_2_2_2xG8_MI300X_160_95p_B1_layer_sp0_patch_size2/checkpoints/paper-STCNNT_HRNET_T1L1G1T1L1G1_T1L1G1T1L1G1_whitney15_hrnet__deployed_T1L1G1T1L1G1_64_ws4_4_16_ps2_2_2_2xG8_MI300X_160_95p_B1_layer_sp0_patch_size2_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_159.pth";     
    "hrnet_100m",  "/data/models/paper_v3/model_search2_hrnet__whitney15_2xG8_MI300X_epoch_160_95p/whitney15_hrnet__deployed_T1L1G1T1L1G1T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_160_95p_B2_layer_sp0/checkpoints/paper-STCNNT_HRNET_T1L1G1T1L1G1T1L1G1T1L1G1_T1L1G1T1L1G1T1L1G1T1L1G1_whitney15_hrnet__deployed_T1L1G1T1L1G1T1L1G1T1L1G1_64_4_4_16_2xG8_MI300X_160_95p_B2_layer_sp0_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_159.pth"; 
    "hrnet_200m",  "/data/models/paper_v3/model_search2_hrnet__whitney15_4xG8_MI300X_epoch_162_98p/whitney15_hrnet__deployed_200m_2dt_4xG8/checkpoints/paper-STCNNT_HRNET_T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1_T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1_whitney15_hrnet__deployed_200m_2dt_4xG8_4x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/best_checkpoint_epoch_137.pth";     

    "hrnet_TTT",                "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1_T1T1T1_kings01_hrnet_2__T1T1T1_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "hrnet_TTTTTT",             "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_T1T1T1T1T1T1_T1T1T1T1T1T1_kings01_hrnet_2__T1T1T1T1T1T1_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "hrnet_Swin",              "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_2xG16_V100_80_90p_B1/checkpoints/paper-STCNNT_HRNET_S3ShS3ShS3Sh_S3ShS3ShS3Sh_kings01_hrnet_2__S3ShS3ShS3Sh_64_8_8_16_2xG16_V100_80_90p_B1_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "hrnet_ViT3D",              "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__V3V3V3_64_8_8_16_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_V3V3V3_V3V3V3_kings01_hrnet_2__V3V3V3_64_8_8_16_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    "hrnet_C3",                 "/data/models/paper_v3/model_search2_hrnet_2__kings01_2xG16_V100_epoch_80_90p/kings01_hrnet_2__C3C3C3_64_2xG16_V100_80_90p_B2/checkpoints/paper-STCNNT_HRNET_C3C3C3_C3C3C3_kings01_hrnet_2__C3C3C3_64_2xG16_V100_80_90p_B2_2x32G16-V100_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_79.pth"; 
    };

model_list = {
    "hrnet_28m", "hrnet_28m.pth";
    "hrnet_50m",  "hrnet_50m.pth";     
    "hrnet_100m", "hrnet_100m.pth"; 
    "hrnet_200m", "hrnet_200m.pth";   

    "hrnet_TTT", "hrnet_TTT.pth";
    "hrnet_TTTTTT", "hrnet_TTTTTT.pth";
    "hrnet_Swin", "hrnet_Swin.pth";
    "hrnet_ViT3D", "hrnet_ViT3D.pth";
    "hrnet_ViT3D_large", "hrnet_ViT3D_large.pth";
    "hrnet_C3", "hrnet_C3.pth";
    "hrnet_C3_large", "hrnet_C3_large.pth"

    %"hrnet_28m_v2", "/data/models/paper_v3/model_search2_hrnet__whitney15_2xG8_MI300X_epoch_160_95p/whitney15_hrnet__deployed_T1L1G1_2xG8/checkpoints/paper2-STCNNT_HRNET_T1L1G1_T1L1G1_whitney15_hrnet__deployed_T1L1G1_2xG8_2x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/best_checkpoint_epoch_156.pth"
    };

% model_list = {
%     "hrnet_200m", "hrnet_200m.pth";   
%     };

model_list_used = model_list;
if ~isempty(input_model_list)
    model_list_used = input_model_list;
end

model_list_used

[fpath, fname, ext] = fileparts(run_file_name);

rec = [];
ind = 1;
for i=1:size(case_lists, 1)
    a_case = case_lists{i, 1};        
    if ~exist(a_case)
        a_case = fullfile(data_dir, a_case);
    end
    for k=1:size(model_list_used, 1)
        a_model_str = model_list_used{k, 1};
        if rep > 1
            a_model_str = strjoin([a_model_str '_rep' num2str(rep)], '');
        end
        a_model = model_list_used{k, 2};
        if ~exist(a_model)
            a_model = fullfile(model_dir, a_model);
        end
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

        a_local_output_case_dir = fullfile(local_output_dir, case_lists{i, 1}, a_model_str);
        a_output_case_dir = fullfile(output_dir, case_lists{i, 1}, a_model_str);

        if check_processed
            if exist(fullfile(a_local_output_case_dir, 'output_imag.npy')) & exist(fullfile(a_local_output_case_dir, 'output_real.npy'))
                continue
            end
            if exist(fullfile(a_local_output_case_dir, 'output.npy'))
                continue
            end
        end

        cmd_str = create_run_command(code_dir, a_case, a_output_case_dir, a_model_str, a_model, no_gmap, no_MR_noise, no_snr_unit, dicom_mode, added_noise_sd, rep, batch_size, check_processed);
        
        rec = [rec; {cmd_str}];
        ind = ind + 1;
    end 
end

num_cases = size(rec, 1);
B = linspace(1, num_cases, nnodes+1)
B = round(B)

for n=1:nnodes

    rec_used = rec(B(n):B(n+1)-1);

    fname_used = [fname '_node_' num2str(n-1)]

    for g=1:numel(gpu_lists)
        if isnumeric(gpu_lists(g))
            gpu_str = num2str(gpu_lists(g));
        else
            gpu_str = gpu_lists(g);
        end
    
        a_run_name = fullfile(fpath, [fname_used '_gpus_' num2str(g) ext ]);
    
        f = fopen(a_run_name, "w")
        fprintf(f, '%s\n', "#!/bin/bash");
        fprintf(f, '%s\n', strjoin(["export CUDA_VISIBLE_DEVICES=", gpu_str], ""));
        fprintf(f, '%s\n', "export DISABLE_FLOAT16_INFERENCE=True");
        fprintf(f, '%s\n', "export RES_DIR=model_res");
        fprintf(f, '%s\n', 'echo "------------------------------------------------------------------------"');
    
        for kk=g:numel(gpu_lists):numel(rec_used)
            cmd_str = rec_used{kk};
            fprintf(f, cmd_str);
            fprintf(f, '\n#############################################\n\n');
        end    
    
        fprintf(f, '\n echo "All run completed ... "\n\n');
        fclose(f);
        
        disp(['run file is - ' a_run_name])
    end
end

end

function cmd_str = create_run_command(code_dir, a_case, a_output_case_dir, a_model_str, a_model, no_gmap, no_MR_noise, no_snr_unit, dicom_mode, added_noise_sd, rep, batch_size, check_processed)

    if rep <= 1
        code_str = strjoin([code_dir "/projects/mri_imaging/inference/run_inference.py --check_processed " num2str(check_processed)], "");
        cmd_str = ["python3 " code_str " --input_dir "  a_case " --output_dir " a_output_case_dir " --scaling_factor 1.0 --im_scaling 1.0 --gmap_scaling 1.0 --input_fname input --gmap_fname gmap --saved_model_path " a_model " --batch_size " batch_size];
    else
        code_str = strjoin([code_dir "/projects/mri_imaging/inference/run_inference_snr_pseudo_replica.py --check_processed " num2str(check_processed)], "");
        cmd_str = ["python3 " code_str " --input_dir "  a_case " --output_dir " a_output_case_dir " --scaling_factor 1.0 --im_scaling 1.0 --gmap_scaling 1.0 --input_fname input --gmap_fname gmap --saved_model_path " a_model " --batch_size " batch_size " --rep " num2str(rep) " --added_noise_sd "  num2str(added_noise_sd)];
    end
    if (no_snr_unit)
        cmd_str = [cmd_str " --scale_by_signal"];
    end
    if (no_gmap)
        cmd_str = [cmd_str " --no_gmap --scaling_factor 1.0"];
    end
    if (dicom_mode)
        cmd_str = [cmd_str " --no_gmap --scale_by_signal"];
    end
    if ~isempty(strfind(a_model_str, "ViT"))
        cmd_str = [cmd_str " --pad_time"];
    end
    if ~isempty(strfind(a_model_str, "Swin"))
        cmd_str = [cmd_str " --pad_time"];
    end

    %cmd_str = [cmd_str " --no_gmap"];
    cmd_str = strjoin(cmd_str, " ")
end