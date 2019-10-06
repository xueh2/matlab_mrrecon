
function seq_params = PerformGadgetronRecon_Matlab_FlowMapping_Get_Seq_Params(dataDir, h5Name, resDir)
% seq_params = PerformGadgetronRecon_Matlab_FlowMapping_Get_Seq_Params(dataDir, h5Name, resDir)

UTDir = getenv('GTPLUS_UT_DIR');

aif_N_runup = 3;
aif_seq_type = 'Flash';
aif_Gd_method = 'LUT';

FA_PD = 5;
N_runup = 3;
seq_type = 'SSFP';
seq_type_PD = 'Flash';
Gd_method = 'LUT';

T1_0_blood_3T = 2000;
T2_0_blood_3T = 250;
T1_0_myo_3T = 1300;
T2_0_myo_3T = 45;

T1_0_blood_1p5T = 1700;
T2_0_blood_1p5T = 220;
T1_0_myo_1p5T = 1100;
T2_0_myo_1p5T = 45;

study_dates = h5Name;
is_h5 = 1;
try
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(h5Name);
    res_dir = fullfile(resDir, study_dates, h5Name)
catch
    res_dir = fullfile(resDir, h5Name)
    is_h5 = 0;
end

HCT = 0.42;
hct_str = num2str(HCT);
ind = find(hct_str=='.');
if(~isempty(ind))
    hct_str(ind(:)) = 'p';
end

%% read in h5 file
if(is_h5)
    dset = ismrmrd.Dataset(fullfile(dataDir, study_dates, [h5Name '.h5']));
else
    dset = ismrmrd.Dataset(fullfile(resDir, h5Name, [h5Name '.h5']));
end
header = ismrmrd.xml.deserialize(dset.readxml());
dset.close();

FA_Perf = header.sequenceParameters.flipAngle_deg(1)
NumOfPD = header.userParameters.userParameterLong(2).value

NumOfConcatenations = 1;

for ii=1:numel(header.userParameters.userParameterLong)
    if(strcmp(header.userParameters.userParameterLong(ii).name, 'NumOfConcatenations')==1)
        NumOfConcatenations = header.userParameters.userParameterLong(ii).value;
        break;
    end
    
    if(strcmp(header.userParameters.userParameterLong(ii).name, 'NumOfProtonDensityImages')==1)
        NumOfPD = header.userParameters.userParameterLong(ii).value;
        break;
    end
end

seq_type = header.sequenceParameters.sequence_type
header.acquisitionSystemInformation.systemFieldStrength_T

if(header.acquisitionSystemInformation.systemFieldStrength_T>2)
    r1 = 4.918800;
    r2 = 5.913300;
    if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
        FA_Perf = FA_Perf * ((header.userParameters.userParameterDouble(6).value / header.userParameters.userParameterDouble(5).value) / 0.8436)
    else
        FA_Perf = FA_Perf * ((header.userParameters.userParameterDouble(6).value / header.userParameters.userParameterDouble(5).value) / 0.44433)
    end

    T1_0_blood = T1_0_blood_3T
    T2_0_blood = T2_0_blood_3T
    T1_0_myo = T1_0_myo_3T
    T2_0_myo = T2_0_myo_3T
else
    r1 = 5.565;
    r2 = 5.78;

    T1_0_blood = T1_0_blood_1p5T
    T2_0_blood = T2_0_blood_1p5T
    T1_0_myo = T1_0_myo_1p5T
    T2_0_myo = T2_0_myo_1p5T
end

if(FA_Perf>header.sequenceParameters.flipAngle_deg(1))
    FA_Perf = header.sequenceParameters.flipAngle_deg(1)
end

aif_echo_time_0 = header.sequenceParameters.TE(2)
aif_echo_time_1 = header.sequenceParameters.TE(3)
aif_FA_PD = header.sequenceParameters.flipAngle_deg(2)
aif_FA_Perf = header.sequenceParameters.flipAngle_deg(2)
if(header.sequenceParameters.TI(2)==0)
    aif_TR = 2.12
else
    aif_TR = header.sequenceParameters.echo_spacing(2)
end
aif_E1_full = header.encoding(2).encodingLimits.kspace_encoding_step_1.maximum + 1
aif_accel_factor = header.encoding(2).parallelImaging.accelerationFactor.kspace_encoding_step_1

if(header.sequenceParameters.TI(2)==0)
    aif_TD = 4.7
    aif_TS = aif_TD + ceil( floor(aif_E1_full / aif_accel_factor) / 2.0 ) * aif_TR
else
    aif_TS = header.sequenceParameters.TI(2)
    aif_TD = aif_TS - ceil( floor(aif_E1_full / aif_accel_factor) / 2.0 ) * aif_TR
end

FA_PD = 5
TR = header.sequenceParameters.echo_spacing(1)
E1_full = header.encoding(1).encodedSpace.matrixSize.y
accel_factor = header.encoding(1).parallelImaging.accelerationFactor.kspace_encoding_step_1

if(header.sequenceParameters.TI(2)==0)
    if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
        TD = 17
        TS = TD + TR * ceil( floor(E1_full / accel_factor) / 2.0) - N_runup * TR
    else
        TD = 8.5
        TS = TD + TR * ceil( floor(E1_full / accel_factor) / 2.0);
    end
else
    TS = header.sequenceParameters.TI(1)
    if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
        TD = TS - TR * ceil( floor(E1_full / accel_factor) / 2.0) - N_runup * TR
    else
        TD = TS - TR * ceil( floor(E1_full / accel_factor) / 2.0)
    end
end

SLC = header.encoding(1).encodingLimits.slice.maximum+1;

%% load images
% 
% cd(fullfile(res_dir, 'DebugOutput'))
% 
% aif_0 = analyze75read('aif_for_TwoEcho_T2StartCorrection_0.hdr');
% aif_cin = analyze75read('aif_cin.hdr');
% aif_cin_Gd = analyze75read('aif_cin_all_echo0_LUTCorrection.hdr');
% aif_cin_Gd_without_R2Star = analyze75read('cin_all_echo0_without_R2Star_LUTCorrection.hdr');
% aif_cin_Gd_baseline_corrected = analyze75read('aif_cin_echo0_all_signal_baseline_corrected.hdr');
% aif_cin_all_echo0_signal = analyze75read('aif_cin_all_echo0_signal.hdr');
% aif_cin_all_echo1_signal = analyze75read('aif_cin_all_echo1_signal.hdr');
% aif_cin_all_echo0_signal_after_R2StarCorrection = analyze75read('aif_cin_all_echo0_signal_after_R2StarCorrection.hdr');
% aif_cin_all_echo0_OverPD_after_R2StarCorrection = analyze75read('aif_cin_all_echo0_OverPD_after_R2StarCorrection.hdr');
% aif_cin_all_R2Star = analyze75read('aif_cin_all_R2Star.hdr');
% aif_cin_all_R2Star_SLEP = analyze75read('aif_cin_all_R2Star_SLEP.hdr');
% cin_used_all_for_computeFlowMap = analyze75read('cin_used_all_for_computeFlowMap.hdr');
% aif_PD = analyze75read('aifPD_for_TwoEcho_T2StartCorrection_0.hdr');
% aif_mask = analyze75read('aif_LV_mask_for_TwoEcho_T2StartCorrection_0.hdr');
% aif_mask_final = analyze75read('AifLVMask_after_Picking.hdr');
% PDMaskIntensity = analyze75read('PDMaskIntensity.hdr');
% pdPicked = analyze75read('pdPicked.hdr');
% aifPicked = analyze75read('aifPicked.hdr');
% % aifMaskIntensity = analyze75read('aifMaskIntensity.hdr');
% % aif0 = analyze75read('AIF_input_for_moco_0_MAG.hdr');
% aif_Gd = analyze75read('aif_cin_Gd_Valid.hdr');
% aif_LUT = analyze75read('aif_cin_LUT_Valid.hdr');
% % aif_PD_echo0 = analyze75read('input_aif_PD_0.hdr');
% % aif_PD_echo1 = analyze75read('input_aif_PD_1.hdr');
% dstAcqTimes = analyze75read('dstAcqTimes_0.hdr');
% AIF_AcqTimes = analyze75read('AIF_AcqTimes_0.hdr');
% perf_LUT = analyze75read('Perf_T1_Correction_LUT.hdr');
% aif_cin_foot_peak_valley = analyze75read('aif_cin_foot_peak_valley.hdr');

seq_params = ws2struct;
