
function [ROITable, sf, rf, sf_i, rf_i, res_table] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, excluded, processing_always, reviewFlag, prefix_SNR, pause_cases, processing_snr_always)
% [ROITable, sf, rf, sf_i, rf_i, res_table] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, excluded, processing_always, reviewFlag, prefix_SNR, pause_cases)
% [ROITable, sf, rf, sf_i, rf_i] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area', 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area_ROI')
% [ROITable, sf, rf, sf_i, rf_i] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area_Native_T1_0', 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area_ROI')
% [ROITable, sf, rf, sf_i, rf_i] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, 'I:\ReconResults\BARTS', 'D:\data\ut\NewData\PaperResults\Barts_ROI')

if(nargin<8)
    excluded = 0;
end

if(nargin<9)
    processing_always = 1;
end

if(nargin<10)
    reviewFlag = 0;
end

if(nargin<11)
    prefix_SNR = 'Perf_Image';
end

if(nargin<12)
    pause_cases = 0;
end

if(nargin<13)
    processing_snr_always = 0;
end

ROITable = PerfTable(1, :);
ROITable = [ROITable 'stress flow' 'rest flow' 'stress ischemia flow' 'rest ischemia flow' 'stress aif peak Gd with baseline correction' 'rest aif peak Gd with baseline correction' 'hematocrit' ...
    'stress E' 'rest E' 'stress PS' 'rest PS' ... 
    'stress Visf' 'rest Visf' 'stress Vp' 'rest Vp' ... 
    'stress Ki_MF' 'rest Ki_MF' 'stress Ki_Fermi' 'rest Ki_Fermi' ... 
    'stress Ki_TwoCompExp' 'rest Ki_TwoCompExp' 'stress Ki_BTEX' 'rest Ki_BTEX' ... 
    'stress E_i' 'rest E_i' 'stress PS_i' 'rest PS_i' ... 
    'stress Visf_i' 'rest Visf_i' 'stress Vp_i' 'rest Vp_i' ... 
    'stress Ki_MF_i' 'rest Ki_MF_i' 'stress Ki_Fermi_i' 'rest Ki_Fermi_i' ... 
    'stress Ki_TwoCompExp_i' 'rest Ki_TwoCompExp_i' 'stress Ki_BTEX_i' 'rest Ki_BTEX_i' ];

nV = numel(PerfTable(1, :));

scanInd = [];
patientID = [];
scanDate = [];
scanTime = [];
age = [];
gender = [];
stressHB = [];
restHB = [];
hematocrit = [];

rest_time = [];
stress_time = [];

sf = [];
rf = [];
sf_i = [];
rf_i = [];

sE = [];
rE = [];
sE_i = [];
rE_i = [];

sPS = [];
rPS = [];
sPS_i = [];
rPS_i = [];

sVisf = [];
rVisf = [];
sVisf_i = [];
rVisf_i = [];

sVp = [];
rVp = [];
sVp_i = [];
rVp_i = [];

sKi_MF = [];
rKi_MF = [];
sKi_MF_i = [];
rKi_MF_i = [];

sKi_Fermi = [];
rKi_Fermi = [];
sKi_Fermi_i = [];
rKi_Fermi_i = [];

sKi_TwoCompExp = [];
rKi_TwoCompExp = [];
sKi_TwoCompExp_i = [];
rKi_TwoCompExp_i = [];

sKi_BTEX = [];
rKi_BTEX = [];
sKi_BTEX_i = [];
rKi_BTEX_i = [];

sDelay = [];
rDelay = [];
sDelay_pixel = [];
rDelay_pixel = [];

sf_sd = [];
sVisf_sd = [];
sVp_sd = [];
sPS_sd = [];
sf_sd_m = [];
sVisf_sd_m = [];
sVp_sd_m = [];
sPS_sd_m = [];

rf_sd = [];
rVisf_sd = [];
rVp_sd = [];
rPS_sd = [];
rf_sd_m = [];
rVisf_sd_m = [];
rVp_sd_m = [];
rPS_sd_m = [];

sf_mean = [];
rf_mean = [];

sE_mean = [];
rE_mean = [];

sVp_mean = [];
rVp_mean = [];

sVisf_mean = [];
rVisf_mean = [];

sPS_mean = [];
rPS_mean = [];

sKi_MF_mean = [];
rKi_MF_mean = [];
sKi_Fermi_mean = [];
rKi_Fermi_mean = [];
sKi_TwoCompExp_mean = [];
rKi_TwoCompExp_mean = [];
sKi_BTEX_mean = [];
rKi_BTEX_mean = [];

sSD_mean = [];
rSD_mean = [];
sPS_SD_mean = [];
rPS_SD_mean = [];
sVisf_SD_mean = [];
rVisf_SD_mean = [];
sVp_SD_mean = [];
rVp_SD_mean = [];
sCC_F_PS_mean = [];
rCC_F_PS_mean = [];
sCC_F_Vp_mean = [];
rCC_F_Vp_mean = [];
sCC_F_Visf_mean = [];
rCC_F_Visf_mean = [];
sCC_PS_Vp_mean = [];
rCC_PS_Vp_mean = [];
sCC_PS_Visf_mean = [];
rCC_PS_Visf_mean = [];
sCC_Vp_Visf_mean = [];
rCC_Vp_Visf_mean = [];
    
sSD = [];
rSD = [];
sSD_i = [];
rSD_i = [];

sPS_SD = [];
rPS_SD = [];
sVisf_SD = [];
rVisf_SD = [];
sVp_SD = [];
rVp_SD = [];

% 4 by 4, F, PS, Vp, Visf
sCC_F_PS = [];
rCC_F_PS = [];
sCC_F_Vp = [];
rCC_F_Vp = [];
sCC_F_Visf = [];
rCC_F_Visf = [];
sCC_PS_Vp = [];
rCC_PS_Vp = [];
sCC_PS_Visf = [];
rCC_PS_Visf = [];
sCC_Vp_Visf = [];
rCC_Vp_Visf = [];

sSNR = [];
rSNR = [];

sGd = [];
rGd = [];

sf_SD = [];
sVisf_SD = [];
sVp_SD = [];
sPS_SD = [];
sE_SD = [];

sf_pixel = [];
sVisf_pixel = [];
sVp_pixel = [];
sPS_pixel = [];
sE_pixel = [];
sKi_MF_pixel = [];
sKi_Fermi_pixel = [];
sKi_TwoCompExp_pixel = [];

sf_delay = [];
sf_delay_pixel = [];

rf_SD = [];
rVisf_SD = [];
rVp_SD = [];
rPS_SD = [];
rE_SD = [];

rf_pixel = [];
rVisf_pixel = [];
rVp_pixel = [];
rPS_pixel = [];
rE_pixel = [];
rKi_MF_pixel = [];
rKi_Fermi_pixel = [];
rKi_TwoCompExp_pixel = [];
rf_delay = [];
rf_delay_pixel = [];

pre_T1_blood = [];
pre_T1_myo = [];

post_T1_blood = [];
post_T1_myo = [];

unused_slices = [];
            
num_column = size(PerfTable, 2);
num = size(PerfTable, 1)-1;

rest_aif_peak = zeros(num,1);
rest_aif_peak_no_T2Star = zeros(num,1);
rest_aif_peak_intensity = zeros(num,1);
rest_aif_peak_intensity_no_T2Star = zeros(num,1);

rest_aif_valley = zeros(num,1);
rest_aif_valley_no_T2Star = zeros(num,1);
rest_aif_valley_intensity = zeros(num,1);
rest_aif_valley_intensity_no_T2Star = zeros(num,1);

rest_aif_T2S_peak = zeros(num, 1);
rest_aif_T2S_baseline = zeros(num, 1);

rest_aif_auc = zeros(num,1);
rest_aif_duration = zeros(num,1);
HeartRate_rest = zeros(num, 1);

stress_aif_valley = zeros(num,1);
stress_aif_valley_no_T2Star = zeros(num,1);
stress_aif_valley_intensity = zeros(num,1);
stress_aif_valley_intensity_no_T2Star = zeros(num,1);

stress_aif_peak = zeros(num,1);
stress_aif_peak_no_T2Star = zeros(num,1);
stress_aif_peak_intensity = zeros(num,1);
stress_aif_peak_intensity_no_T2Star = zeros(num,1);

stress_aif_T2S_peak = zeros(num, 1);
stress_aif_T2S_baseline = zeros(num, 1);

stress_aif_auc = zeros(num,1);
stress_aif_duration = zeros(num,1);
HeartRate_stress = zeros(num, 1);

ecv = [];

for n=1:num
       
    ii = n;
    
    stressCase = PerfTable{n+1, stress_column};
    restCase = PerfTable{n+1, rest_column};
    if(ischemia_column>0 & ischemia_column<=num_column)
        Category = PerfTable{n+1, ischemia_column};
    else
        Category = 1;
    end
    
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]); 
    disp(['==================================================================']);  
       
    [configName, scannerID, pID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(stressCase);
    [configName, scannerID, pID, studyID, measurementID, study_dates_rest, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(restCase);

    figDir = fullfile(resDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_' study_time_stress '_' study_time_rest '_Figure']) 
    
    if(excluded)
        roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_ROI_excluded'])
    else
        roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_ROI'])
    end
    
    SNRResult_file = fullfile(roiDir, [prefix_SNR '_SNR_PerfResult_' scannerID '_' pID '_' studyID '_' study_dates '.mat']);
    Gd_Result_file = fullfile(roiDir, [prefix_SNR '_Gd_PerfResult_' scannerID '_' pID '_' studyID '_' study_dates '.mat']);

    stressDir = fullfile(resDir, study_dates, stressCase)
    restDir = fullfile(resDir, study_dates, restCase)
    
    scanInd = [scanInd; n];
    patientID = [patientID; {pID}];
    scanDate = [scanDate; study_dates];
    scanTime = [scanTime; study_time_stress];
    if(size(PerfTable, 2)>9)
        stressHB = [stressHB; PerfTable{n+1, 9}];
        restHB = [restHB; PerfTable{n+1, 13}];
        age = [age; PerfTable{n+1, 14}];
        gender = [gender; PerfTable{n+1, 15}];
    else
        stressHB = [stressHB; -1];
        restHB = [restHB; -1];
        age = [age; PerfTable{n+1, 6}];
        gender = [gender; PerfTable{n+1, 5}];
    end
    
    if(hct_column>0)
        hematocrit = [hematocrit; PerfTable{n+1, hct_column}];

        HCT = PerfTable{n+1, hct_column};
        if(isnan(HCT))
            HCT = 0;
        else
            if(HCT>1)
                HCT = HCT/100;
            end
        end
    else
        hematocrit = [hematocrit; 0.42];
        HCT = 0.42;
    end
    
    %% aif related
    try
        [HeartRate_rest(n), aif_cin_Gd_rest, aif_cin_Gd_rest_without_R2Star, aif_cin_all_echo0_signal_rest, aif_cin_all_echo0_signal_after_R2StarCorrection_rest, footTime_rest, peakTime_rest, valleyTime_rest, auc_rest, R2Star_rest, SampledInterval_rest, KiMap, flowMap, EMap, PSMap, VisfMap, VpMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, restCase);
        has_rest = 1;
    catch
        aif_cin_Gd_rest = -1;
        aif_cin_Gd_rest_without_R2Star = -1;
        aif_cin_all_echo0_signal_rest = -1;
        aif_cin_all_echo0_signal_after_R2StarCorrection_rest= -1; 
        footTime_rest = -1;
        peakTime_rest = -1;
        valleyTime_rest = -1;
        R2Star_rest = -1;
        SampledInterval_rest = -1;
        has_rest = 0;
    end
    
    [HeartRate_stress(n), aif_cin_Gd_stress, aif_cin_Gd_stress_without_R2Star, aif_cin_all_echo0_signal_stress, aif_cin_all_echo0_signal_after_R2StarCorrection_stress, footTime_stress, peakTime_stress, valleyTime_stress, auc_stress, R2Star_stress, SampledInterval_stress, KiMap, flowMap, EMap, PSMap, VisfMap, VpMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, stressCase);

    rest_time = [rest_time; {study_time_rest}];
    stress_time = [stress_time; {study_time_stress}];

    rest_aif_peak(n) = max(aif_cin_Gd_rest);
    rest_aif_peak_no_T2Star(n) = max(aif_cin_Gd_rest_without_R2Star);
    rest_aif_peak_intensity(n) = max(aif_cin_all_echo0_signal_after_R2StarCorrection_rest);
    rest_aif_peak_intensity_no_T2Star(n) = max(aif_cin_all_echo0_signal_rest);

    if(has_rest)
        r1 = ceil(valleyTime_rest);
        r2 = ceil(valleyTime_rest-0.5);
        r3 = ceil(valleyTime_rest+0.5);

        rest_aif_valley(n) = (aif_cin_Gd_rest(r1) + aif_cin_Gd_rest(r2) + aif_cin_Gd_rest(r3))/3;
        rest_aif_valley_no_T2Star(n) = (aif_cin_Gd_rest_without_R2Star(r1) + aif_cin_Gd_rest_without_R2Star(r2) + aif_cin_Gd_rest_without_R2Star(r3))/3;
        rest_aif_valley_intensity(n) = (aif_cin_all_echo0_signal_after_R2StarCorrection_rest(r1) + aif_cin_all_echo0_signal_after_R2StarCorrection_rest(r2) + aif_cin_all_echo0_signal_after_R2StarCorrection_rest(r3))/3;
        rest_aif_valley_intensity_no_T2Star(n) = (aif_cin_all_echo0_signal_rest(r1) + aif_cin_all_echo0_signal_rest(r2) + aif_cin_all_echo0_signal_rest(r3))/3;

        rest_aif_auc(n) = auc_rest;
        
        rest_aif_duration(n) = (valleyTime_rest - footTime_rest) * 0.5;

        figure; hold on; plot(aif_cin_all_echo0_signal_rest); plot(aif_cin_all_echo0_signal_after_R2StarCorrection_rest, 'r'); hold off
        figure; hold on; plot(aif_cin_Gd_rest); plot(aif_cin_Gd_rest_without_R2Star, 'r');

        rest_T2S = 1.0 ./ (R2Star_rest+eps);
        stress_T2S = 1.0 ./ (R2Star_stress+eps);

        rest_aif_T2S_peak(n) = min( [rest_T2S(ceil(peakTime_rest)) rest_T2S(ceil(peakTime_rest+0.5)) rest_T2S(ceil(peakTime_rest-0.5))] );
        stress_aif_T2S_peak(n) = min( [stress_T2S(ceil(peakTime_stress)) stress_T2S(ceil(peakTime_stress+0.5)) stress_T2S(ceil(peakTime_stress-0.5))] );

        v = rest_T2S(3:floor(footTime_rest));
        ind = find(v<100);
        if (~isempty(ind)) rest_aif_T2S_baseline(n) = mean(v(ind)); end

        v = stress_T2S(3:floor(footTime_stress));
        ind = find(v<100);
        if (~isempty(ind)) stress_aif_T2S_baseline(n) = mean(v(ind)); end
    end
    % -----------------------------------

    stress_aif_peak(n) = max(aif_cin_Gd_stress);
    stress_aif_peak_no_T2Star(n) = max(aif_cin_Gd_stress_without_R2Star);
    stress_aif_peak_intensity(n) = max(aif_cin_all_echo0_signal_after_R2StarCorrection_stress);
    stress_aif_peak_intensity_no_T2Star(n) = max(aif_cin_all_echo0_signal_stress);

    r1 = ceil(valleyTime_stress);
    r2 = ceil(valleyTime_stress-0.5);
    r3 = ceil(valleyTime_stress+0.5);

    stress_aif_valley(n) = (aif_cin_Gd_stress(r1) + aif_cin_Gd_stress(r2) + aif_cin_Gd_stress(r3))/3;
    stress_aif_valley_no_T2Star(n) = (aif_cin_Gd_stress_without_R2Star(r1) + aif_cin_Gd_stress_without_R2Star(r2) + aif_cin_Gd_stress_without_R2Star(r3))/3;
    stress_aif_valley_intensity(n) = (aif_cin_all_echo0_signal_after_R2StarCorrection_stress(r1) + aif_cin_all_echo0_signal_after_R2StarCorrection_stress(r2) + aif_cin_all_echo0_signal_after_R2StarCorrection_stress(r3))/3;
    stress_aif_valley_intensity_no_T2Star(n) = (aif_cin_all_echo0_signal_stress(r1) + aif_cin_all_echo0_signal_stress(r2) + aif_cin_all_echo0_signal_stress(r3))/3;

    stress_aif_auc(n) = auc_stress;

    stress_aif_duration(n) = (valleyTime_stress - footTime_stress) * SampledInterval_rest/1e3;
            
    closeall
    %% mapping related
    
    if(isFileExist(fullfile(roiDir, 'myo_stress1.mat')) | isFileExist(fullfile(roiDir, 'myo_stress2.mat')) | isFileExist(fullfile(roiDir, 'myo_stress3.mat')))
        s1_roi = 'myo_stress1.mat';
        s2_roi = 'myo_stress2.mat';
        s3_roi = 'myo_stress3.mat';
        
        r1_roi = 'myo_rest1.mat';
        r2_roi = 'myo_rest2.mat';
        r3_roi = 'myo_rest3.mat';
    else
        s1_roi = 's1.mat';
        s2_roi = 's2.mat';
        s3_roi = 's3.mat';
        
        r1_roi = 'r1.mat';
        r2_roi = 'r2.mat';
        r3_roi = 'r3.mat';        
    end
        
    s1_visf_roi = 's_visf_1.mat';
    s2_visf_roi = 's_visf_2.mat';
    s3_visf_roi = 's_visf_3.mat';

    r1_visf_roi = 'r_visf_1.mat';
    r2_visf_roi = 'r_visf_2.mat';
    r3_visf_roi = 'r_visf_3.mat';        

    has_roi_1 = isFileExist(fullfile(roiDir, s1_roi)) | isFileExist(fullfile(roiDir, r1_roi));
    has_roi_2 = isFileExist(fullfile(roiDir, s2_roi)) | isFileExist(fullfile(roiDir, r2_roi));
    has_roi_3 = isFileExist(fullfile(roiDir, s3_roi)) | isFileExist(fullfile(roiDir, r3_roi));
    
    if(has_roi_1 | has_roi_2 | has_roi_3)
       disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]);       
        
       two_ROI = 0;

        PerfResult_file = fullfile(figDir, 'PerfResult.mat');
        if(~processing_always & isFileExist(PerfResult_file))
            tt = load(PerfResult_file);
            v = tt.PerfResult;
        else        
            cd(roiDir)
            v = PerfTable(n+1, :);

            nV = numel(v);

            if(isFileExist(fullfile(roiDir, s1_roi)))
                s1 = load(fullfile(roiDir, s1_roi));
            else
                s1 = [];
            end
            
            if(isFileExist(fullfile(roiDir, s2_roi)))
                s2 = load(fullfile(roiDir, s2_roi));
            else
                s2 = [];
            end
            
            
            if(isFileExist(fullfile(roiDir, s3_roi)))
                s3 = load(fullfile(roiDir, s3_roi));
            else
                s3 = [];
            end

            if(isFileExist(fullfile(roiDir, r1_roi)))
                r1 = load(fullfile(roiDir, r1_roi));
            else
                r1 = [];
            end 
            
            if(isFileExist(fullfile(roiDir, r2_roi)))
                r2 = load(fullfile(roiDir, r2_roi));
            else
                r2 = [];
            end
            
            if(isFileExist(fullfile(roiDir, r3_roi)))
                r3 = load(fullfile(roiDir, r3_roi));
            else
                r3 = [];
            end
                       
            if(isFileExist(fullfile(roiDir, s1_visf_roi)))
                s1_visf = load(fullfile(roiDir, s1_visf_roi));
            else
                s1_visf = [];
            end
            
            if(isFileExist(fullfile(roiDir, s2_visf_roi)))
                s2_visf = load(fullfile(roiDir, s2_visf_roi));
            else
                s2_visf = [];
            end
            
            
            if(isFileExist(fullfile(roiDir, s3_visf_roi)))
                s3_visf = load(fullfile(roiDir, s3_visf_roi));
            else
                s3_visf = [];
            end

            if(isFileExist(fullfile(roiDir, r1_visf_roi)))
                r1_visf = load(fullfile(roiDir, r1_visf_roi));
            else
                r1_visf = [];
            end 
            
            if(isFileExist(fullfile(roiDir, r2_visf_roi)))
                r2_visf = load(fullfile(roiDir, r2_visf_roi));
            else
                r2_visf = [];
            end
            
            if(isFileExist(fullfile(roiDir, r3_visf_roi)))
                r3_visf = load(fullfile(roiDir, r3_visf_roi));
            else
                r3_visf = [];
            end
            
            if(has_rest)
                rest = load(fullfile(figDir, 'rest.mat'));
                res_rest = PerformGadgetronRecon_SavedIsmrmrd_ROIValues_OneCase(r1, r2, r3, rest, r1_visf, r2_visf, r3_visf);
            else
                res_rest.flow = [-1 -1 -1];
                res_rest.E = [-1 -1 -1];
                res_rest.Visf = [-1 -1 -1];
                res_rest.Vp = [-1 -1 -1];
                res_rest.PS = [-1 -1 -1];
                res_rest.SD = [-1 -1 -1];
                res_rest.E_i = [-1 -1 -1];
                res_rest.Ki_MF = [-1 -1 -1];
                res_rest.Ki_Fermi = [-1 -1 -1];
                res_rest.Ki_TwoCompExp = [-1 -1 -1];
                res_rest.Ki_BTEX = [-1 -1 -1];
                rest.aif_rest_baseline_corrected = -1;
                res_rest.delay = [-1 -1 -1];
            end
            
            stress = load(fullfile(figDir, 'stress.mat'));
            res_stress = PerformGadgetronRecon_SavedIsmrmrd_ROIValues_OneCase(s1, s2, s3, stress, s1_visf, s2_visf, s3_visf);
            
            v{nV+1} = res_stress.flow;
            disp(['stress flow : ' num2str(v{nV+1})]);
            disp(['stress Ki_Fermi : ' num2str(res_stress.Ki_Fermi)]);
            
            v{nV+2} = res_rest.flow;
            disp(['rest flow : ' num2str(v{nV+2})]);

            has_rest_second_roi = 0;
            if( (~isempty(s1) & numel(s1.ROI_info_table)==2) | (~isempty(s2) & numel(s2.ROI_info_table)==2) | (~isempty(s3) & numel(s3.ROI_info_table)==2))
                
                two_ROI = 1;

                v{nV+3} = res_stress.flow_i;
                disp(['stress flow, ischemia : ' num2str(v{nV+3})]);
                
                if(isfield(res_rest, 'flow_i'))
                    v{nV+4} = res_rest.flow_i;
                    has_rest_second_roi = 1;
                else
                    v{nV+4} = [-1 -1 -1];
                end
            else
                v{nV+3} = [-1 -1 -1];
                v{nV+4} = [-1 -1 -1];
            end
            
            v{nV+5} = max(stress.aif_stress_baseline_corrected);
            v{nV+6} = max(rest.aif_rest_baseline_corrected);
            v{nV+7} = HCT;
            
            
            if(reviewFlag)
                
                try
                    if(~isempty(s1)) figure; imagescn(stress.flow_stress(:,:,1,end), [0 8], [], [], [], fullfile(roiDir, s1_roi)); PerfColorMap; end
                    if(~isempty(s2)) figure; imagescn(stress.flow_stress(:,:,2,end), [0 8], [], [], [], fullfile(roiDir, s2_roi)); PerfColorMap; end
                    if(~isempty(s3)) figure; imagescn(stress.flow_stress(:,:,3,end), [0 8], [], [], [], fullfile(roiDir, s3_roi)); PerfColorMap; end

                    if(~isempty(s1)) figure; imagescn(stress.Vp_stress(:,:,1,end), [0 20], [], [], [], fullfile(roiDir, s1_roi)); MBVColorMap; end
                    if(~isempty(s2)) figure; imagescn(stress.Vp_stress(:,:,2,end), [0 20], [], [], [], fullfile(roiDir, s2_roi)); MBVColorMap; end
                    if(~isempty(s3)) figure; imagescn(stress.Vp_stress(:,:,3,end), [0 20], [], [], [], fullfile(roiDir, s3_roi)); MBVColorMap; end

                    if(~isempty(s1_visf)) figure; imagescn(stress.Visf_stress(:,:,1,end), [0 80], [], [], [], fullfile(roiDir, s1_visf_roi)); ECVColorMap; end
                    if(~isempty(s2_visf)) figure; imagescn(stress.Visf_stress(:,:,2,end), [0 80], [], [], [], fullfile(roiDir, s2_visf_roi)); ECVColorMap; end
                    if(~isempty(s3_visf)) figure; imagescn(stress.Visf_stress(:,:,3,end), [0 80], [], [], [], fullfile(roiDir, s3_visf_roi)); ECVColorMap; end

                    if(has_rest)
                        if(~isempty(r1)) figure; imagescn(rest.flow_rest(:,:,1,end), [0 8], [], [], [], fullfile(roiDir, r1_roi)); PerfColorMap; end
                        if(~isempty(r2)) figure; imagescn(rest.flow_rest(:,:,2,end), [0 8], [], [], [], fullfile(roiDir, r2_roi)); PerfColorMap; end
                        if(~isempty(r3)) figure; imagescn(rest.flow_rest(:,:,3,end), [0 8], [], [], [], fullfile(roiDir, r3_roi)); PerfColorMap; end

                        if(~isempty(r1)) figure; imagescn(rest.Vp_rest(:,:,1,end), [0 20], [], [], [], fullfile(roiDir, r1_roi)); MBVColorMap; end
                        if(~isempty(r2)) figure; imagescn(rest.Vp_rest(:,:,2,end), [0 20], [], [], [], fullfile(roiDir, r2_roi)); MBVColorMap; end
                        if(~isempty(r3)) figure; imagescn(rest.Vp_rest(:,:,3,end), [0 20], [], [], [], fullfile(roiDir, r3_roi)); MBVColorMap; end

                        if(~isempty(r1_visf)) figure; imagescn(rest.Visf_rest(:,:,1,end), [0 80], [], [], [], fullfile(roiDir, r1_visf_roi)); ECVColorMap; end
                        if(~isempty(r2_visf)) figure; imagescn(rest.Visf_rest(:,:,2,end), [0 80], [], [], [], fullfile(roiDir, r2_visf_roi)); ECVColorMap; end
                        if(~isempty(r3_visf)) figure; imagescn(rest.Visf_rest(:,:,3,end), [0 80], [], [], [], fullfile(roiDir, r3_visf_roi)); ECVColorMap; end
                    end
                catch
                end
                
                % plot the histogram
                try
                    figure;

                    subplot(4, 2, 1)
                    pt = Perfusion_GetResultValusPerPixel(res_stress.flow_pixels);
                    ind = find(pt>0 & pt<5.5);
                    histfit(pt(ind), 64, 'lognormal');
                    xlabel('Stress flow, ml/min/g')
                    xlim([0 8])
                    subplot(4, 2, 2)
                    pt = Perfusion_GetResultValusPerPixel(res_rest.flow_pixels);
                    ind = find(pt>0);
                    histfit(pt(ind), 64, 'lognormal');
                    xlabel('Rest flow, ml/min/g')
                    xlim([0 8])

                    subplot(4, 2, 3)
                    pt = Perfusion_GetResultValusPerPixel(res_stress.PS_pixels);
                    [N, X] = hist(pt, 64);
                    ind = find(pt>0.42);
                    histfit(pt(ind), 64, 'lognormal');
                    xlabel('Stress PS, ml/min/g')
                    xlim([0 8])
                    subplot(4, 2, 4)
                    pt = Perfusion_GetResultValusPerPixel(res_rest.PS_pixels);
                    ind = find(pt>0);
                    histfit(pt(ind), 64, 'lognormal');
                    xlabel('Rest PS, ml/min/g')
                    xlim([0 8])

                    subplot(4, 2, 5)
                    pt = Perfusion_GetResultValusPerPixel(res_stress.Visf_pixels);
                    ind = find(pt>0);
                    histfit(pt(ind), 64, 'lognormal');
                    xlabel('Stress Visf, ml/g')
                    xlim([0 70])
                    subplot(4, 2, 6)
                    pt = Perfusion_GetResultValusPerPixel(res_rest.Visf_pixels);
                    ind = find(pt>0);
                    histfit(pt(ind), 64, 'lognormal');
                    xlabel('Rest Visf, ml/g')
                    xlim([0 70])

                    subplot(4, 2, 7)
                    pt = Perfusion_GetResultValusPerPixel(res_stress.Vp_pixels);
                    ind = find(pt>0);
                    histfit(pt(ind), 64, 'kernel');
                    xlabel('Stress Vb, ml/g')
                    xlim([0 20])
                    subplot(4, 2, 8)
                    pt = Perfusion_GetResultValusPerPixel(res_rest.Vp_pixels);
                    ind = find(pt>0);
                    histfit(pt(ind), 64, 'kernel');
                    xlabel('Rest Vb, ml/g')
                    xlim([0 20])
                catch
                end
                
                if(pause_cases) 
                    user_in = input('unaccepted slice for stress:', 's');
                    if(numel(user_in)>0)
                        unaccepted_slice = str2num(user_in);
                        save('unused_slices_stress.mat', 'unaccepted_slice');
                    end
                    
                    user_in = input('unaccepted slice for rest:', 's');
                    if(numel(user_in)>0)
                        unaccepted_slice = str2num(user_in);
                        save('unused_slices_rest.mat', 'unaccepted_slice');
                    end
                    
                    user_in = input('accept cases y or n :');
                    if(user_in=='n')
                        onlyReview = 1;
                        [h_flow_stress, h_flow_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir,  stressCase, restCase, [0 6], onlyReview, resDir);      
                        pause;
                    end
                end
                closeall
            end

            if(isFileExist('unused_slices_stress.mat'))
                u_slc = load('unused_slices_stress.mat');                
            else
                u_slc.unaccepted_slice = [];
            end
            
            if(isFileExist('unused_slices_rest.mat'))
                u_slc2 = load('unused_slices_rest.mat');                
            else
                u_slc2.unaccepted_slice = [];
            end
            
            unused_slices =[unused_slices; {u_slc.unaccepted_slice} {u_slc2.unaccepted_slice}];
            
            ind = 8;
            
            % E
            v{nV+ind} = res_stress.E; ind = ind+1;
            v{nV+ind} = res_rest.E; ind = ind+1;
            
            % PS
            v{nV+ind} = res_stress.PS; ind = ind+1;
            v{nV+ind} = res_rest.PS; ind = ind+1;
            
            % Visf
            v{nV+ind} = res_stress.Visf; ind = ind+1;
            v{nV+ind} = res_rest.Visf; ind = ind+1;
            
            % Vp
            v{nV+ind} = res_stress.Vp; ind = ind+1;
            v{nV+ind} = res_rest.Vp; ind = ind+1;
            
            % Ki_MF
            v{nV+ind} = res_stress.Ki_MF; ind = ind+1;
            v{nV+ind} = res_rest.Ki_MF; ind = ind+1;
                        
            % Ki_Fermi
            v{nV+ind} = res_stress.Ki_Fermi; ind = ind+1;
            v{nV+ind} = res_rest.Ki_Fermi; ind = ind+1;
            
            % Ki_TwoCompExp
            v{nV+ind} = res_stress.Ki_TwoCompExp; ind = ind+1;
            v{nV+ind} = res_rest.Ki_TwoCompExp; ind = ind+1;
            
            % Ki_BTEX
            v{nV+ind} = res_stress.Ki_BTEX; ind = ind+1;
            v{nV+ind} = res_rest.Ki_BTEX; ind = ind+1;
            
            sE = [sE; res_stress.E];
            rE = [rE; res_rest.E];

            sPS = [sPS; res_stress.PS];
            rPS = [rPS; res_rest.PS];

            sVisf = [sVisf; res_stress.Visf];
            rVisf = [rVisf; res_rest.Visf];

            sVp = [sVp; res_stress.Vp];
            rVp = [rVp; res_rest.Vp];

            sKi_MF = [sKi_MF; res_stress.Ki_MF];
            rKi_MF = [rKi_MF; res_rest.Ki_MF];

            sKi_Fermi = [sKi_Fermi; res_stress.Ki_Fermi];
            rKi_Fermi = [rKi_Fermi; res_rest.Ki_Fermi];

            sKi_TwoCompExp = [sKi_TwoCompExp; res_stress.Ki_TwoCompExp];
            rKi_TwoCompExp = [rKi_TwoCompExp; res_rest.Ki_TwoCompExp];

            sKi_BTEX = [sKi_BTEX; res_stress.Ki_BTEX];
            rKi_BTEX = [rKi_BTEX; res_rest.Ki_BTEX];

            sSD = [sSD; res_stress.SD];
            rSD = [rSD; res_rest.SD];
            sPS_SD = [sPS_SD; res_stress.PS_SD];
            rPS_SD = [rPS_SD; res_rest.PS_SD];
            sVisf_SD = [sVisf_SD; res_stress.Visf_SD];
            rVisf_SD = [rVisf_SD; res_rest.Visf_SD];
            sVp_SD = [sVp_SD; res_stress.Vp_SD];
            rVp_SD = [rVp_SD; res_rest.Vp_SD];
            
            sCC_F_PS = [sCC_F_PS; res_stress.CC_F_PS];
            rCC_F_PS = [rCC_F_PS; res_rest.CC_F_PS];
            
            sCC_F_Vp = [sCC_F_Vp; res_stress.CC_F_Vp];
            rCC_F_Vp = [rCC_F_Vp; res_rest.CC_F_Vp];
            
            sCC_F_Visf = [sCC_F_Visf; res_stress.CC_F_Visf];
            rCC_F_Visf = [rCC_F_Visf; res_rest.CC_F_Visf];
            
            sCC_PS_Vp = [sCC_PS_Vp; res_stress.CC_PS_Vp];
            rCC_PS_Vp = [rCC_PS_Vp; res_rest.CC_PS_Vp];
            
            sCC_PS_Visf = [sCC_PS_Visf; res_stress.CC_PS_Visf];
            rCC_PS_Visf = [rCC_PS_Visf; res_rest.CC_PS_Visf];
            
            sCC_Vp_Visf = [sCC_Vp_Visf; res_stress.CC_Vp_Visf];
            rCC_Vp_Visf = [rCC_Vp_Visf; res_rest.CC_Vp_Visf];
            
            sDelay = [sDelay; res_stress.delay];
            rDelay = [rDelay; res_rest.delay];
            sDelay_pixel = [sDelay_pixel; res_stress.delay_pixels];
            rDelay_pixel = [rDelay_pixel; res_rest.delay_pixels];

            sf_pixel        = [sf_pixel; res_stress.flow_pixels];
            sE_pixel        = [sE_pixel; res_stress.E_pixels];
            sVisf_pixel     = [sVisf_pixel; res_stress.Visf_pixels];
            sVp_pixel       = [sVp_pixel; res_stress.Vp_pixels];
            sPS_pixel       = [sPS_pixel; res_stress.PS_pixels];
            sKi_MF_pixel    = [sKi_MF_pixel; res_stress.Ki_MF_pixels];
            sKi_Fermi_pixel = [sKi_Fermi_pixel; res_stress.Ki_Fermi_pixels];
            sKi_TwoCompExp_pixel = [sKi_TwoCompExp_pixel; res_stress.Ki_TwoCompExp_pixels];
            
            sf_sd = [sf_sd; res_stress.flow_sd];
            sPS_sd = [sPS_sd; res_stress.PS_sd];
            sVisf_sd = [sVisf_sd; res_stress.Visf_sd];
            sVp_sd = [sVp_sd; res_stress.Vp_sd];
            
            sf_delay = [sf_delay; {res_stress.delay_flow}];
            sf_delay_pixel = [sf_delay_pixel; {res_stress.delay_flow_pixels}];
            
            rf_pixel        = [rf_pixel; res_rest.flow_pixels];
            rE_pixel        = [rE_pixel; res_rest.E_pixels];
            rVisf_pixel     = [rVisf_pixel; res_rest.Visf_pixels];
            rVp_pixel       = [rVp_pixel; res_rest.Vp_pixels];
            rPS_pixel       = [rPS_pixel; res_rest.PS_pixels];
            rKi_MF_pixel    = [rKi_MF_pixel; res_rest.Ki_MF_pixels];
            rKi_Fermi_pixel = [rKi_Fermi_pixel; res_rest.Ki_Fermi_pixels];
            rKi_TwoCompExp_pixel = [rKi_TwoCompExp_pixel; res_rest.Ki_TwoCompExp_pixels];
                        
            rf_sd = [rf_sd; res_rest.flow_sd];
            rPS_sd = [rPS_sd; res_rest.PS_sd];
            rVisf_sd = [rVisf_sd; res_rest.Visf_sd];
            rVp_sd = [rVp_sd; res_rest.Vp_sd];
            
            rf_delay = [rf_delay; {res_rest.delay_flow}];
            rf_delay_pixel = [rf_delay_pixel; {res_rest.delay_flow_pixels}];

            if(two_ROI)
                % E
                v{nV+ind} = res_stress.E_i; ind = ind+1;
                v{nV+ind} = res_rest.E_i; ind = ind+1;

                % PS
                v{nV+ind} = res_stress.PS_i; ind = ind+1;
                v{nV+ind} = res_rest.PS_i; ind = ind+1;

                % Visf
                v{nV+ind} = res_stress.Visf_i; ind = ind+1;
                v{nV+ind} = res_rest.Visf_i; ind = ind+1;

                % Vp
                v{nV+ind} = res_stress.Vp_i; ind = ind+1;
                v{nV+ind} = res_rest.Vp_i; ind = ind+1;

                % Ki_MF
                v{nV+ind} = res_stress.Ki_MF_i; ind = ind+1;
                v{nV+ind} = res_rest.Ki_MF_i; ind = ind+1;

                % Ki_Fermi
                v{nV+ind} = res_stress.Ki_Fermi_i; ind = ind+1;
                v{nV+ind} = res_rest.Ki_Fermi_i; ind = ind+1;

                % Ki_TwoCompExp
                v{nV+ind} = res_stress.Ki_TwoCompExp_i; ind = ind+1;
                v{nV+ind} = res_rest.Ki_TwoCompExp_i; ind = ind+1;

                % Ki_BTEX
                v{nV+ind} = res_stress.Ki_BTEX_i; ind = ind+1;
                v{nV+ind} = res_rest.Ki_BTEX_i; ind = ind+1;
                
                sE_i = [sE_i; res_stress.E_i];
                rE_i = [rE_i; res_rest.E_i];

                sPS_i = [sPS_i; res_stress.PS_i];
                rPS_i = [rPS_i; res_rest.PS_i];

                sVisf_i = [sVisf_i; res_stress.Visf_i];
                rVisf_i = [rVisf_i; res_rest.Visf_i];

                sVp_i = [sVp_i; res_stress.Vp_i];
                rVp_i = [rVp_i; res_rest.Vp_i];

                sKi_MF_i = [sKi_MF_i; res_stress.Ki_MF_i];
                rKi_MF_i = [rKi_MF_i; res_rest.Ki_MF_i];

                sKi_Fermi_i = [sKi_Fermi_i; res_stress.Ki_Fermi_i];
                rKi_Fermi_i = [rKi_Fermi_i; res_rest.Ki_Fermi_i];

                sKi_TwoCompExp_i = [sKi_TwoCompExp_i; res_stress.Ki_TwoCompExp_i];
                rKi_TwoCompExp_i = [rKi_TwoCompExp_i; res_rest.Ki_TwoCompExp_i];

                sKi_BTEX_i = [sKi_BTEX_i; res_stress.Ki_BTEX_i];
                rKi_BTEX_i = [rKi_BTEX_i; res_rest.Ki_BTEX_i];

                sSD_i = [sSD_i; res_stress.SD_i];
                rSD_i = [rSD_i; res_rest.SD_i];
            else
                % E
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % PS
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Visf
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Vp
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Ki_MF
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Ki_Fermi
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Ki_TwoCompExp
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Ki_BTEX
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                
                sE_i = [sE_i; -1 -1 -1];
                rE_i = [rE_i; -1 -1 -1];

                sPS_i = [sPS_i; -1 -1 -1];
                rPS_i = [rPS_i; -1 -1 -1];

                sVisf_i = [sVisf_i; -1 -1 -1];
                rVisf_i = [rVisf_i; -1 -1 -1];

                sVp_i = [sVp_i; -1 -1 -1];
                rVp_i = [rVp_i; -1 -1 -1];

                sKi_MF_i = [sKi_MF_i; -1 -1 -1];
                rKi_MF_i = [rKi_MF_i; -1 -1 -1];

                sKi_Fermi_i = [sKi_Fermi_i; -1 -1 -1];
                rKi_Fermi_i = [rKi_Fermi_i; -1 -1 -1];

                sKi_TwoCompExp_i = [sKi_TwoCompExp_i; -1 -1 -1];
                rKi_TwoCompExp_i = [rKi_TwoCompExp_i; -1 -1 -1];

                sKi_BTEX_i = [sKi_BTEX_i; -1 -1 -1];
                rKi_BTEX_i = [rKi_BTEX_i; -1 -1 -1];

                sSD_i = [sSD_i; -1 -1 -1];
                rSD_i = [rSD_i; -1 -1 -1];
            end
            
            % ---------------------------------------------
            % load T1 map
            t1_pre = fullfile(roiDir, 'T1_pre.mat');
            t1_post = fullfile(roiDir, 'T1_post.mat');
            
            pre_t1_blood = [0 0 0];
            pre_t1_myo = [0 0 0];
            post_t1_blood = [0 0 0];
            post_t1_myo = [0 0 0];
            
            if(isFileExist(t1_pre))
                t1_v = load(t1_pre);
                pre_t1_blood = [t1_v.ROI_info_table(1,1).ROI_mean t1_v.ROI_info_table(3,2).ROI_mean t1_v.ROI_info_table(5,3).ROI_mean];
                pre_t1_myo = [t1_v.ROI_info_table(2,1).ROI_mean t1_v.ROI_info_table(4,2).ROI_mean t1_v.ROI_info_table(6,3).ROI_mean];
            end
            
            if(isFileExist(t1_post))
                t1_v = load(t1_post);
                post_t1_blood = [t1_v.ROI_info_table(1,1).ROI_mean t1_v.ROI_info_table(3,2).ROI_mean t1_v.ROI_info_table(5,3).ROI_mean];
                post_t1_myo = [t1_v.ROI_info_table(2,1).ROI_mean t1_v.ROI_info_table(4,2).ROI_mean t1_v.ROI_info_table(6,3).ROI_mean];
            end
            
            disp(['pre_t1_blood is ' num2str(pre_t1_blood)])
            disp(['post_t1_blood is ' num2str(post_t1_blood)])
            disp(['pre_t1_myo is ' num2str(pre_t1_myo)])
            disp(['post_t1_myo is ' num2str(post_t1_myo)])

            pre_T1_blood = [pre_T1_blood; pre_t1_blood];
            pre_T1_myo = [pre_T1_myo; pre_t1_myo];
            
            post_T1_blood = [post_T1_blood; post_t1_blood];
            post_T1_myo = [post_T1_myo; post_t1_myo];
            
            ecv = [ecv; (1-HCT) * (1./post_t1_myo - 1./pre_t1_myo) ./ (1./post_t1_blood - 1./pre_t1_blood)];
            
            % ---------------------------------------------
            % SNR
            
            % if(processing_snr_always | ~isFileExist(SNRResult_file))
            if(processing_snr_always)
                try
                    cd(stressDir)
                    s_gfactor = readGTPlusExportImageSeries_Squeeze(300);
                    s_gfactor = flipdim(s_gfactor, 2);

                    moco_perf = readGTPlusExportImageSeries_Squeeze(104);
                    moco_perf = flipdim(moco_perf, 2);
                    
                    sdata1 = squeeze(moco_perf(:,:,1,:));
                    sdata2 = squeeze(moco_perf(:,:,2,:));
                    sdata3 = squeeze(moco_perf(:,:,3,:));

                    cd(restDir)
                    r_gfactor = readGTPlusExportImageSeries_Squeeze(300);
                    r_gfactor = flipdim(r_gfactor, 2);

                    moco_perf = readGTPlusExportImageSeries_Squeeze(104);
                    moco_perf = flipdim(moco_perf, 2);
                    
                    rdata1 = squeeze(moco_perf(:,:,1,:));
                    rdata2 = squeeze(moco_perf(:,:,2,:));
                    rdata3 = squeeze(moco_perf(:,:,3,:));
                    
                    snr_s1 = 25 * sdata1 ./ squeeze(s_gfactor(:,:,1,:));
                    snr_s2 = 25 * sdata2 ./ squeeze(s_gfactor(:,:,2,:));
                    snr_s3 = 25 * sdata3 ./ squeeze(s_gfactor(:,:,3,:));

                    snr_r1 = 25 * rdata1 ./ squeeze(r_gfactor(:,:,1,:));
                    snr_r2 = 25 * rdata2 ./ squeeze(r_gfactor(:,:,2,:));
                    snr_r3 = 25 * rdata3 ./ squeeze(r_gfactor(:,:,3,:));

                    cd(roiDir)
                    if(~isempty(s1))
                        BW1=zeros(size(sdata1(:,:,1)));
                        BW1=roipoly(sdata1(:,:,1), s1(1).ROI_info_table(1).ROI_x_coordinates, s1(1).ROI_info_table(1).ROI_y_coordinates);
                        index1=find(BW1 >0);
                    end

                    if(~isempty(s2))
                        BW2=zeros(size(sdata2(:,:,1)));
                        BW2=roipoly(sdata2(:,:,1), s2(1).ROI_info_table(1).ROI_x_coordinates, s2(1).ROI_info_table(1).ROI_y_coordinates);
                        index2=find(BW2 >0);
                    end

                    if(~isempty(s3))
                        BW3=zeros(size(sdata3(:,:,1)));
                        BW3=roipoly(sdata3(:,:,1), s3(1).ROI_info_table(1).ROI_x_coordinates, s3(1).ROI_info_table(1).ROI_y_coordinates);
                        index3=find(BW3 >0);
                    end

                    if(~isempty(r1))
                        rBW1=zeros(size(rdata1(:,:,1)));
                        rBW1=roipoly(rdata1(:,:,1), r1(1).ROI_info_table(1).ROI_x_coordinates, r1(1).ROI_info_table(1).ROI_y_coordinates);
                        rindex1=find(rBW1 >0);
                    end

                    if(~isempty(r2))
                        rBW2=zeros(size(rdata2(:,:,1)));
                        rBW2=roipoly(rdata2(:,:,1), r2(1).ROI_info_table(1).ROI_x_coordinates, r2(1).ROI_info_table(1).ROI_y_coordinates);
                        rindex2=find(rBW2 >0);
                    end

                    if(~isempty(r3))
                        rBW3=zeros(size(rdata3(:,:,1)));
                        rBW3=roipoly(rdata3(:,:,1), r3(1).ROI_info_table(1).ROI_x_coordinates, r3(1).ROI_info_table(1).ROI_y_coordinates);
                        rindex3=find(rBW3 >0);
                    end

                    % gfactor scaled by 100
                    % data is scaled by 4
                    nRep = size(sdata1, 3);

                    s_SNR = zeros(nRep, 3);
                    r_SNR = zeros(nRep, 3);

                    for rr=1:nRep
                        sd1 = snr_s1(:,:,rr);
                        sd2 = snr_s2(:,:,rr);
                        sd3 = snr_s3(:,:,rr);

                        rd1 = snr_r1(:,:,rr);
                        rd2 = snr_r2(:,:,rr);
                        rd3 = snr_r3(:,:,rr);

                        v1 = 0;
                        if(~isempty(s1))
                            v1 = mean(sd1(index1));
                        end

                        v2 = 0;
                        if(~isempty(s2))
                            v2 = mean(sd2(index2));
                        end

                        v3 = 0;
                        if(~isempty(s3))
                            v3 = mean(sd3(index3));
                        end

                        s_SNR(rr, :) = [v1 v2 v3];

                        v1 = 0;
                        if(~isempty(r1))
                            v1 = mean(rd1(rindex1));
                        end

                        v2 = 0;
                        if(~isempty(r2))
                            v2 = mean(rd2(rindex2));
                        end

                        v3 = 0;
                        if(~isempty(r3))
                            v3 = mean(rd3(rindex3));
                        end
                        r_SNR(rr, :) = [v1 v2 v3];                
                    end

                    sSNR = [sSNR; max(s_SNR(4:end, :))];
                    rSNR = [rSNR; max(r_SNR(4:end, :))];

                    s_gfactor_map = s_gfactor(:,:,:,1) /100;
                    r_gfactor_map = r_gfactor(:,:,:,1) /100;
                    save(SNRResult_file, 's_SNR', 'r_SNR', 'snr_s1', 'snr_s2', 'snr_s3', 'snr_r1', 'snr_r2', 'snr_r3', 's_gfactor_map', 'r_gfactor_map');               
                catch
                    s_SNR = 0;
                    r_SNR = 0;
                    sSNR = [sSNR; -1 -1 -1];
                    rSNR = [rSNR; -1 -1 -1];
                end
            else
                try
                snr_v = load(SNRResult_file);
                
                s_SNR = snr_v.s_SNR;
                r_SNR = snr_v.r_SNR;
                
                sSNR = [sSNR; max(s_SNR(4:end, :))];
                rSNR = [rSNR; max(r_SNR(4:end, :))];
                catch
                    sSNR = [sSNR; -1 -1 -1];
                    rSNR = [rSNR; -1 -1 -1];
                end
            end
        end
    
        ROITable = [ROITable; v];
                
        rf = [rf; ROITable{end, nV+2}];         
        
        if(two_ROI)
            sf_i = [sf_i; ROITable{end, nV+3}];
            rf_i = [rf_i; ROITable{end, nV+4}];
        else
            sf_i = [sf_i; -1 -1 -1];
            rf_i = [rf_i; -1 -1 -1];
        end
        sf = [sf; ROITable{end, nV+1}];
        
        PerfResult = v;
        save(fullfile(figDir, 'PerfResult.mat'), 'PerfResult');
        
        closeall
    else
        disp(['Cannot find the ROI ... ']);
        winopen(figDir);
        mkdir(roiDir);
        winopen(roiDir);
        % pause;
        
        sf = [sf; -1 -1 -1];
        rf = [rf; -1 -1 -1];
        sf_i = [sf_i; -1 -1 -1];
        rf_i = [rf_i; -1 -1 -1];

        sE = [sE; -1 -1 -1];
        rE = [rE; -1 -1 -1];
        sE_i = [sE_i; -1 -1 -1];
        rE_i = [rE_i; -1 -1 -1];

        sPS = [sPS; -1 -1 -1];
        rPS = [rPS; -1 -1 -1];
        sPS_i = [sPS_i; -1 -1 -1];
        rPS_i = [rPS_i; -1 -1 -1];

        sVisf = [sVisf; -1 -1 -1];
        rVisf = [rVisf; -1 -1 -1];
        sVisf_i = [sVisf_i; -1 -1 -1];
        rVisf_i = [rVisf_i; -1 -1 -1];

        sVp = [sVp; -1 -1 -1];
        rVp = [rVp; -1 -1 -1];
        sVp_i = [sVp_i; -1 -1 -1];
        rVp_i = [rVp_i; -1 -1 -1];

        sKi_MF = [sKi_MF; -1 -1 -1];
        rKi_MF = [rKi_MF; -1 -1 -1];
        sKi_MF_i = [sKi_MF_i; -1 -1 -1];
        rKi_MF_i = [rKi_MF_i; -1 -1 -1];

        sKi_Fermi = [sKi_Fermi; -1 -1 -1];
        rKi_Fermi = [rKi_Fermi; -1 -1 -1];
        sKi_Fermi_i = [sKi_Fermi_i; -1 -1 -1];
        rKi_Fermi_i = [rKi_Fermi_i; -1 -1 -1];

        sKi_TwoCompExp = [sKi_TwoCompExp; -1 -1 -1];
        rKi_TwoCompExp = [rKi_TwoCompExp; -1 -1 -1];
        sKi_TwoCompExp_i = [sKi_TwoCompExp_i; -1 -1 -1];
        rKi_TwoCompExp_i = [rKi_TwoCompExp_i; -1 -1 -1];

        sKi_BTEX = [sKi_BTEX; -1 -1 -1];
        rKi_BTEX = [rKi_BTEX; -1 -1 -1];
        sKi_BTEX_i = [sKi_BTEX_i; -1 -1 -1];
        rKi_BTEX_i = [rKi_BTEX_i; -1 -1 -1];

        sDelay = [sDelay; -1 -1 -1];
        rDelay = [rDelay; -1 -1 -1];

%         sf_mean = [sf_mean; -1];
%         rf_mean = [rf_mean; -1];
% 
%         sE_mean = [sE_mean; -1];
%         rE_mean = [rE_mean; -1];
% 
%         sVp_mean = [sVp_mean; -1];
%         rVp_mean = [rVp_mean; -1];
% 
%         sVisf_mean = [sVisf_mean; -1];
%         rVisf_mean = [rVisf_mean; -1];
% 
%         sPS_mean = [sPS_mean; -1];
%         rPS_mean = [rPS_mean; -1];
% 
%         sKi_MF_mean = [sKi_MF_mean; -1];
%         rKi_MF_mean = [rKi_MF_mean; -1];
%         sKi_Fermi_mean = [sKi_Fermi_mean; -1];
%         rKi_Fermi_mean = [rKi_Fermi_mean; -1];
%         sKi_TwoCompExp_mean = [sKi_TwoCompExp_mean; -1];
%         rKi_TwoCompExp_mean = [rKi_TwoCompExp_mean; -1];
%         sKi_BTEX_mean = [sKi_BTEX_mean; -1];
%         rKi_BTEX_mean = [rKi_BTEX_mean; -1];

        sSD = [sSD; -1 -1 -1];
        rSD = [rSD; -1 -1 -1];
        sSD_i = [sSD_i; -1 -1 -1];
        rSD_i = [rSD_i; -1 -1 -1];

        sPS_SD = [sPS_SD; -1 -1 -1];
        rPS_SD = [rPS_SD; -1 -1 -1];
        sVisf_SD = [sVisf_SD; -1 -1 -1];
        rVisf_SD = [rVisf_SD; -1 -1 -1];
        sVp_SD = [sVp_SD; -1 -1 -1];
        rVp_SD = [rVp_SD; -1 -1 -1];

        sCC_F_PS = [sCC_F_PS; -1 -1 -1];
        rCC_F_PS = [sCC_F_PS; -1 -1 -1];

        sCC_F_Vp = [sCC_F_Vp; -1 -1 -1];
        rCC_F_Vp = [rCC_F_Vp; -1 -1 -1];

        sCC_F_Visf = [sCC_F_Visf; -1 -1 -1];
        rCC_F_Visf = [rCC_F_Visf; -1 -1 -1];

        sCC_PS_Vp = [sCC_PS_Vp; -1 -1 -1];
        rCC_PS_Vp = [rCC_PS_Vp; -1 -1 -1];

        sCC_PS_Visf = [sCC_PS_Visf; -1 -1 -1];
        rCC_PS_Visf = [rCC_PS_Visf; -1 -1 -1];

        sCC_Vp_Visf = [sCC_Vp_Visf; -1 -1 -1];
        rCC_Vp_Visf = [rCC_Vp_Visf; -1 -1 -1];
            
        sSNR = [sSNR; -1 -1 -1];
        rSNR = [rSNR; -1 -1 -1];

        sf_pixel = [sf_pixel; cell(1,3)];
        sVisf_pixel = [sVisf_pixel; cell(1,3)];
        sVp_pixel = [sVp_pixel; cell(1,3)];
        sPS_pixel = [sPS_pixel; cell(1,3)];
        sE_pixel = [sE_pixel; cell(1,3)];
        sKi_MF_pixel = [sKi_MF_pixel; cell(1,3)];

        rf_pixel = [rf_pixel; cell(1,3)];
        rVisf_pixel = [rVisf_pixel; cell(1,3)];
        rVp_pixel = [rVp_pixel; cell(1,3)];
        rPS_pixel = [rPS_pixel; cell(1,3)];
        rE_pixel = [rE_pixel; cell(1,3)];
        rKi_MF_pixel = [rKi_MF_pixel; cell(1,3)];

        pre_T1_blood = [pre_T1_blood; -1 -1 -1];
        pre_T1_myo = [pre_T1_myo; -1 -1 -1];

        post_T1_blood = [post_T1_blood; -1 -1 -1];
        post_T1_myo = [post_T1_myo; -1 -1 -1];
        
        ecv = [ecv; -1 -1 -1];
    end
end

cd(resDir)
save Leeds_res

for n=1:size(sf,1)
    sf_mean = [sf_mean; get_entry_mean(sf(n,:))];
    rf_mean = [rf_mean; get_entry_mean(rf(n,:))];

    sE_mean = [sE_mean; get_entry_mean(sE(n,:))];
    rE_mean = [rE_mean; get_entry_mean(rE(n,:))];
    
    sPS_mean = [sPS_mean; get_entry_mean(sPS(n,:))];
    rPS_mean = [rPS_mean; get_entry_mean(rPS(n,:))];
    
    sVp_mean = [sVp_mean; get_entry_mean(sVp(n,:))];
    rVp_mean = [rVp_mean; get_entry_mean(rVp(n,:))];
    
    sVisf_mean = [sVisf_mean; get_entry_mean(sVisf(n,:))];
    rVisf_mean = [rVisf_mean; get_entry_mean(rVisf(n,:))];

    sKi_MF_mean = [sKi_MF_mean; get_entry_mean(sKi_MF(n,:))];
    rKi_MF_mean = [rKi_MF_mean; get_entry_mean(rKi_MF(n,:))];
    
    sKi_Fermi_mean = [sKi_Fermi_mean; get_entry_mean(sKi_Fermi(n,:))];
    rKi_Fermi_mean = [rKi_Fermi_mean; get_entry_mean(rKi_Fermi(n,:))];
    
    sKi_TwoCompExp_mean = [sKi_TwoCompExp_mean; get_entry_mean(sKi_TwoCompExp(n,:))];
    rKi_TwoCompExp_mean = [rKi_TwoCompExp_mean; get_entry_mean(rKi_TwoCompExp(n,:))];
    
    sKi_BTEX_mean = [sKi_BTEX_mean; get_entry_mean(sKi_BTEX(n,:))];
    rKi_BTEX_mean = [rKi_BTEX_mean; get_entry_mean(rKi_BTEX(n,:))];  
    
    sSD_mean = [sSD_mean; get_entry_mean(sSD(n,:))];
    rSD_mean = [rSD_mean; get_entry_mean(rSD(n,:))];
    sPS_SD_mean = [sPS_SD_mean; get_entry_mean(sPS_SD(n,:))];
    rPS_SD_mean = [rPS_SD_mean; get_entry_mean(rPS_SD(n,:))];
    sVisf_SD_mean = [sVisf_SD_mean; get_entry_mean(sVisf_SD(n,:))];
    rVisf_SD_mean = [rVisf_SD_mean; get_entry_mean(rVisf_SD(n,:))];
    sVp_SD_mean = [sVp_SD_mean; get_entry_mean(sVp_SD(n,:))];
    rVp_SD_mean = [rVp_SD_mean; get_entry_mean(rVp_SD(n,:))];
    sCC_F_PS_mean = [sCC_F_PS_mean; get_entry_mean(sCC_F_PS(n,:))];
    rCC_F_PS_mean = [rCC_F_PS_mean; get_entry_mean(rCC_F_PS(n,:))];
    sCC_F_Vp_mean = [sCC_F_Vp_mean; get_entry_mean(sCC_F_Vp(n,:))];
    rCC_F_Vp_mean = [rCC_F_Vp_mean; get_entry_mean(rCC_F_Vp(n,:))];
    sCC_F_Visf_mean = [sCC_F_Visf_mean; get_entry_mean(sCC_F_Visf(n,:))];
    rCC_F_Visf_mean = [rCC_F_Visf_mean; get_entry_mean(rCC_F_Visf(n,:))];
    sCC_PS_Vp_mean = [sCC_PS_Vp_mean; get_entry_mean(sCC_PS_Vp(n,:))];
    rCC_PS_Vp_mean = [rCC_PS_Vp_mean; get_entry_mean(rCC_PS_Vp(n,:))];
    sCC_PS_Visf_mean = [sCC_PS_Visf_mean; get_entry_mean(sCC_PS_Visf(n,:))];
    rCC_PS_Visf_mean = [rCC_PS_Visf_mean; get_entry_mean(rCC_PS_Visf(n,:))];
    sCC_Vp_Visf_mean = [sCC_Vp_Visf_mean; get_entry_mean(sCC_Vp_Visf(n,:))];
    rCC_Vp_Visf_mean = [rCC_Vp_Visf_mean; get_entry_mean(rCC_Vp_Visf(n,:))];
    
    sf_sd_m = [sf_sd_m; get_entry_mean(sf_sd(n,:))];
    sPS_sd_m = [sPS_sd_m; get_entry_mean(sPS_sd(n,:))];
    sVisf_sd_m = [sVisf_sd_m; get_entry_mean(sVisf_sd(n,:))];
    sVp_sd_m = [sVp_sd_m; get_entry_mean(sVp_sd(n,:))];
    
    rf_sd_m = [rf_sd_m; get_entry_mean(rf_sd(n,:))];
    rPS_sd_m = [rPS_sd_m; get_entry_mean(rPS_sd(n,:))];
    rVisf_sd_m = [rVisf_sd_m; get_entry_mean(rVisf_sd(n,:))];
    rVp_sd_m = [rVp_sd_m; get_entry_mean(rVp_sd(n,:))];
end

sVp_mean
rVp_mean

N = size(scanInd,1);

if(size(age,1)<N)
    age = zeros(N,1);
end

if(size(gender,1)<N)
    gender = zeros(N,1);
end

if(size(hematocrit,1)<N)
    hematocrit = zeros(N,1);
end

res_table = table(scanInd, patientID, scanDate, scanTime, age, gender, stressHB, restHB, hematocrit, ... 
                rest_time, stress_time, ...                
                HeartRate_stress, ... 
                stress_aif_peak, stress_aif_peak_no_T2Star, stress_aif_peak_intensity, stress_aif_peak_intensity_no_T2Star, ... 
                stress_aif_valley, stress_aif_valley_no_T2Star, stress_aif_valley_intensity, stress_aif_valley_intensity_no_T2Star, ... 
                stress_aif_T2S_peak, stress_aif_T2S_baseline, stress_aif_duration, stress_aif_auc, ...                
                HeartRate_rest, ... 
                rest_aif_peak, rest_aif_peak_no_T2Star, rest_aif_peak_intensity, rest_aif_peak_intensity_no_T2Star, ... 
                rest_aif_valley, rest_aif_valley_no_T2Star, rest_aif_valley_intensity, rest_aif_valley_intensity_no_T2Star, ... 
                rest_aif_T2S_peak, rest_aif_T2S_baseline, rest_aif_duration, rest_aif_auc, ...                 
                sf, rf, sf_i, rf_i, ... 
                sE, rE, sE_i, rE_i, ... 
                sPS, rPS, sPS_i, rPS_i, ...
                sVisf, rVisf, sVisf_i, rVisf_i, ... 
                sVp, rVp, sVp_i, rVp_i, ... 
                sKi_MF, rKi_MF, sKi_MF_i, rKi_MF_i, ... 
                sKi_Fermi, rKi_Fermi, sKi_Fermi_i, rKi_Fermi_i, ... 
                sKi_TwoCompExp, rKi_TwoCompExp, sKi_TwoCompExp_i, rKi_TwoCompExp_i, ...
                sKi_BTEX, rKi_BTEX, sKi_BTEX_i, rKi_BTEX_i, ... 
                sSD, rSD, sSD_i, rSD_i, ...
                sPS_SD, rPS_SD, ...
                sVisf_SD, rVisf_SD, ...
                sVp_SD, rVp_SD, ...
                sCC_F_PS, rCC_F_PS, ...
                sCC_F_Vp, rCC_F_Vp, ... 
                sCC_F_Visf, rCC_F_Visf, ...
                sCC_PS_Vp, rCC_PS_Vp, ... 
                sCC_PS_Visf, rCC_PS_Visf, ...
                sCC_Vp_Visf, rCC_Vp_Visf, ...
                sf_mean, rf_mean, ... 
                sE_mean, rE_mean, ... 
                sPS_mean, rPS_mean, ...
                sVisf_mean, rVisf_mean, ... 
                sVp_mean, rVp_mean, ... 
                sKi_MF_mean, rKi_MF_mean, ... 
                sKi_Fermi_mean, rKi_Fermi_mean, ... 
                sKi_TwoCompExp_mean, rKi_TwoCompExp_mean, ...
                sKi_BTEX_mean, rKi_BTEX_mean, ...
                sSD_mean, ...
                rSD_mean, ...
                sPS_SD_mean, ...
                rPS_SD_mean, ...
                sVisf_SD_mean, ...
                rVisf_SD_mean, ...
                sVp_SD_mean, ...
                rVp_SD_mean, ...
                sCC_F_PS_mean, ...
                rCC_F_PS_mean, ...
                sCC_F_Vp_mean, ...
                rCC_F_Vp_mean, ...
                sCC_F_Visf_mean, ...
                rCC_F_Visf_mean, ...
                sCC_PS_Vp_mean, ...
                rCC_PS_Vp_mean, ...
                sCC_PS_Visf_mean, ...
                rCC_PS_Visf_mean, ...
                sCC_Vp_Visf_mean, ...
                rCC_Vp_Visf_mean, ...
                sDelay, rDelay, sDelay_pixel, rDelay_pixel, ..., 
                sf_sd, sf_sd_m, ...
                sPS_sd, sPS_sd_m, ...
                sVisf_sd, sVisf_sd_m, ...
                sVp_sd, sVp_sd_m, ...
                rf_sd, rf_sd_m, ...
                rPS_sd, rPS_sd_m, ...
                rVisf_sd, rVisf_sd_m, ...
                rVp_sd, rVp_sd_m, ...
                sSNR, rSNR, ...
                sf_pixel, sE_pixel, sVisf_pixel, sVp_pixel, sPS_pixel, sKi_MF_pixel, sKi_Fermi_pixel, sKi_TwoCompExp_pixel, ...
                sf_delay, sf_delay_pixel, ...
                rf_pixel, rE_pixel, rVisf_pixel, rVp_pixel, rPS_pixel, rKi_MF_pixel, rKi_Fermi_pixel, rKi_TwoCompExp_pixel, ...                
                rf_delay, rf_delay_pixel, ...
                pre_T1_blood, pre_T1_myo, ...            
                post_T1_blood, post_T1_myo, ...
                ecv, ...
                unused_slices ...
);

disp('=======================================================================');
disp(['Stress flow - ' num2str(mean(sf_mean(:))) '+/-' num2str(std(sf_mean(:)))]);
disp(['Rest flow - ' num2str(mean(rf_mean(:))) '+/-' num2str(std(rf_mean(:)))]);

ind = find(sf_i(:)>0);
disp(['Stress flow, ischemia - ' num2str(mean(sf_i(ind))) '+/-' num2str(std(sf_i(ind)))]);

end

function v = get_entry_mean(f)
    ind = find(f>0);
    if(isempty(ind))
        v=-1;
    else
        v = mean(f(ind));
    end
end
