
function [HeartRate, aif_cin_Gd, aif_cin_Gd_without_R2Star_baseline_corrected, footTime, peakTime, valleyTime, KiMap, flowMap, EMap, PSMap, VisfMap, VpMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, caseName)
% [aif_cin_Gd, aif_cin_Gd_without_R2Star, foot, peak, valley, KiMap, flowMap, EMap, PSMap, VisfMap, BloodVolumeMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, caseName)
% [aif_cin_Gd, aif_cin_Gd_without_R2Star, foot, peak, valley, KiMap, flowMap, EMap, PSMap, VisfMap, BloodVolumeMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan('I:\ReconResults\BARTS', 'Perfusion_AIF_TwoEchoes_Interleaved_R2_42110_196106578_196106587_845_20160613-114338')


[configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(caseName);

debugDir = fullfile(resDir, study_dates, caseName, 'DebugOutput');

aif_0 = analyze75read(fullfile(debugDir, 'aif_for_TwoEcho_T2StartCorrection_0.hdr'));
aif_cin = analyze75read(fullfile(debugDir, 'aif_cin.hdr'));
aif_cin_Gd = analyze75read(fullfile(debugDir, 'aif_cin_all_echo0_LUTCorrection.hdr'));
aif_cin_Gd_without_R2Star = analyze75read(fullfile(debugDir, 'cin_all_echo0_without_R2Star_LUTCorrection.hdr'));
aif_cin_Gd_baseline_corrected = analyze75read(fullfile(debugDir, 'aif_cin_echo0_all_signal_baseline_corrected.hdr'));
aif_cin_all_echo0_signal = analyze75read(fullfile(debugDir, 'aif_cin_all_echo0_signal.hdr'));
aif_cin_all_echo1_signal = analyze75read(fullfile(debugDir, 'aif_cin_all_echo1_signal.hdr'));
aif_cin_all_echo0_signal_after_R2StarCorrection = analyze75read(fullfile(debugDir, 'aif_cin_all_echo0_signal_after_R2StarCorrection.hdr'));
aif_cin_all_echo0_OverPD_after_R2StarCorrection = analyze75read(fullfile(debugDir, 'aif_cin_all_echo0_OverPD_after_R2StarCorrection.hdr'));
aif_cin_all_R2Star = analyze75read(fullfile(debugDir, 'aif_cin_all_R2Star.hdr'));
aif_cin_all_R2Star_SLEP = analyze75read(fullfile(debugDir, 'aif_cin_all_R2Star_SLEP.hdr'));
aif_PD = analyze75read(fullfile(debugDir, 'aifPD_for_TwoEcho_T2StartCorrection_0.hdr'));
aif_mask = analyze75read(fullfile(debugDir, 'aif_LV_mask_for_TwoEcho_T2StartCorrection_0.hdr'));
aif_mask_final = analyze75read(fullfile(debugDir, 'AifLVMask_after_Picking.hdr'));
AIF_AcqTimes = analyze75read(fullfile(debugDir, 'AIF_AcqTimes_0.hdr'));

RRInterval = AIF_AcqTimes(2:end) - AIF_AcqTimes(1:end-1);
HeartRate = 60*1e3/median(RRInterval);
        
sampleinterval = 0.5;
sigmas = [0.8 2.0 3.2];
sigmaMeasure = 0.2;
thresGradient = 0.5;
plotFlag = 1;

[slope_rest, timeToPeak_rest, peakTime, areaUnderCurve_rest, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(aif_cin_Gd_baseline_corrected(:), sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');

aif_rest2 = aif_cin_Gd_baseline_corrected(round(peakTime):end);
[slope, timeToPeak, peakTime2, areaUnderCurve, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(-aif_rest2(:), sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, 0, 'Feature detection of AIF LV signal');   
valleyTime = round(peakTime2) + peakTime - 0.5;
footTime = round(peakTime - timeToPeak_rest);   

v = mean(aif_cin_Gd_without_R2Star(4:footTime));
aif_cin_Gd_without_R2Star_baseline_corrected = aif_cin_Gd_without_R2Star - v;

KiMap = [];
flowMap = [];
EMap = [];
PSMap = [];
VisfMap = [];
VpMap = [];

try
    r1 = analyze75read(fullfile(debugDir, 'flow_maps_after_hole_filling_0.hdr'));
    r2 = analyze75read(fullfile(debugDir, 'flow_maps_after_hole_filling_1.hdr'));
    r3 = analyze75read(fullfile(debugDir, 'flow_maps_after_hole_filling_2.hdr'));
    flow = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
    flow = flipdim(flipdim(flow, 2), 1);
    flowMap = flow;
catch
end

try
    r1 = analyze75read(fullfile(debugDir, 'Ki_maps_after_hole_filling_0.hdr'));
    r2 = analyze75read(fullfile(debugDir, 'Ki_maps_after_hole_filling_1.hdr'));
    r3 = analyze75read(fullfile(debugDir, 'Ki_maps_after_hole_filling_2.hdr'));
    Ki = cat(4, r1, r2, r3);
    Ki = flipdim(flipdim(Ki, 2), 1);
    KiMap = Ki;
catch
end

try
    r1 = analyze75read(fullfile(debugDir, 'PS_maps_after_hole_filling_0.hdr'));
    r2 = analyze75read(fullfile(debugDir, 'PS_maps_after_hole_filling_1.hdr'));
    r3 = analyze75read(fullfile(debugDir, 'PS_maps_after_hole_filling_2.hdr'));
    PS = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
    PS = flipdim(flipdim(PS, 2), 1);
    PSMap = PS;
catch
end

try
    r1 = analyze75read(fullfile(debugDir, 'blood_volume_maps_after_hole_filling_0.hdr'));
    r2 = analyze75read(fullfile(debugDir, 'blood_volume_maps_after_hole_filling_1.hdr'));
    r3 = analyze75read(fullfile(debugDir, 'blood_volume_maps_after_hole_filling_2.hdr'));
    Vp = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
    Vp = flipdim(flipdim(Vp, 2), 1);
    VpMap = Vp;
catch
end

try
    r1 = analyze75read(fullfile(debugDir, 'interstitial_volume_maps_0.hdr'));
    r2 = analyze75read(fullfile(debugDir, 'interstitial_volume_maps_1.hdr'));
    r3 = analyze75read(fullfile(debugDir, 'interstitial_volume_maps_2.hdr'));
    Visf = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
    Visf = flipdim(flipdim(Visf, 2), 1);
    VisfMap = Visf;
catch
end

try
    r1 = analyze75read(fullfile(debugDir, 'E_maps_0.hdr'));
    r2 = analyze75read(fullfile(debugDir, 'E_maps_1.hdr'));
    r3 = analyze75read(fullfile(debugDir, 'E_maps_2.hdr'));
    E = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
    E = flipdim(flipdim(E, 2), 1);
    EMap = E;
catch
end
