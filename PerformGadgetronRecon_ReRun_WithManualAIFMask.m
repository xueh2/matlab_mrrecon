
function PerformGadgetronRecon_ReRun_WithManualAIFMask(perf_case, dataDir, resDir, host, configNamePreset)
% PerformGadgetronRecon_ReRun_WithManualAIFMask(perf_case, dataDir, resDir, host, configNamePreset)

perf_xml = find_perf_xml(configNamePreset, perf_case);
disp([perf_case ' - ' perf_xml]);

[configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(perf_case);

resDir_short = resDir;

resDir = fullfile(resDir, study_dates, perf_case);
cd(resDir)

debug_dir = fullfile(resDir, 'DebugOutput');

aif_cin = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin.hdr'));
aif_cin_Gd = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
aif_cin_Gd_without_R2Star = analyze75read(fullfile(resDir, 'DebugOutput', 'cin_all_echo0_without_R2Star_LUTCorrection.hdr'));
aif_cin_Gd_baseline_corrected = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_echo0_all_signal_baseline_corrected.hdr'));
aif_cin_all_echo0_signal = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_signal.hdr'));
aif_cin_all_echo1_signal = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo1_signal.hdr'));
aif_cin_all_echo0_signal_after_R2StarCorrection = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_signal_after_R2StarCorrection.hdr'));
aif_cin_all_echo0_OverPD_after_R2StarCorrection = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_OverPD_after_R2StarCorrection.hdr'));
aif_cin_all_R2Star = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_R2Star.hdr'));
aif_cin_all_R2Star_SLEP = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_R2Star_SLEP.hdr'));
aif_PD = analyze75read(fullfile(resDir, 'DebugOutput', 'aifPD_for_TwoEcho_T2StartCorrection_0.hdr'));
aif_mask = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_LV_mask_for_TwoEcho_T2StartCorrection_0.hdr'));
aif_mask_final = analyze75read(fullfile(resDir, 'DebugOutput', 'AifLVMask_after_Picking.hdr'));
aif_moco = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_moco.hdr'));

try
    aif_LV_mask_plot = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_LV_mask_plot_.hdr'));
catch
    aif_LV_mask_plot = []
end

aif = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));

figure; imagescn(aif_moco, [], [], [5], 3);
figure; imagescn(aif_mask, [], [], 5);

pause
[mask_file,path] = uigetfile;

roi = load (fullfile(path, mask_file));

BW=zeros(size(aif_moco(:,:,1)));
roi_no = 1;
image_no = 1;
BW=roipoly(aif_moco(:,:,20), roi.ROI_info_table(roi_no, image_no).ROI_x_coordinates, roi.ROI_info_table(roi_no, image_no).ROI_y_coordinates);
figure; imagescn(BW, [], [], 5);

index=find(BW >0);
BW(index) = 1;

header = CreateGtImageHeader(BW);
Matlab_gt_write_analyze(single(BW), header, fullfile(path, 'manual_mask'));

mkdir('d:/temp/aif_masks')
copyfile(fullfile(path, 'manual_mask.hdr'), fullfile('d:/temp/aif_masks', [perf_case '_manual_mask.hdr']));
copyfile(fullfile(path, 'manual_mask.img'), fullfile('d:/temp/aif_masks', [perf_case '_manual_mask.img']));

config_dir = getenv('GT_CONDIG_DIR');
cd(config_dir)

if(strcmp(host, 'localhost')==0)
    cd('d:/temp/aif_masks')
    command = ['scp ./' perf_case '_manual_mask.hdr xueh2@' host ':/tmp/gadgetron_data/'];
    dos(command)
    command = ['scp ./' perf_case '_manual_mask.img xueh2@' host ':/tmp/gadgetron_data/'];
    dos(command)
    
    cd(config_dir)

    command = ['sed -i "s%AIF_LV_MaskFile</name><value></value>%' 'AIF_LV_MaskFile</name><value>/tmp/gadgetron_data/' perf_case '_manual_mask</value>' '%g" ' perf_xml]
    dos(command);
else
    command = ['sed -i "s%AIF_LV_MaskFile</name><value></value>%' 'AIF_LV_MaskFile</name><value>d:/temp/aif_masks/' perf_case '_manual_mask</value>' '%g" ' perf_xml]
    dos(command);
end

checkProcessed = 0;
delete_old_res = 1;
startRemoteGT = 1;

PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, {perf_case}, host, resDir_short, checkProcessed, delete_old_res, startRemoteGT, {perf_xml});

if(strcmp(host, 'localhost')==0)
    cd(config_dir)

    command = ['sed -i "s%AIF_LV_MaskFile</name><value>/tmp/gadgetron_data/' perf_case '_manual_mask</value>%' 'AIF_LV_MaskFile</name><value></value>' '%g" ' perf_xml]
    dos(command);
else
    command = ['sed -i "s%AIF_LV_MaskFile</name><value>d:/temp/aif_masks/' perf_case '_manual_mask</value>%' 'AIF_LV_MaskFile</name><value></value>' '%g" ' perf_xml]
    dos(command);
end
