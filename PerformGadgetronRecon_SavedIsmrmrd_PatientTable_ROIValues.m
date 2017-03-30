
function [ROITable, sf, rf, sf_i, rf_i, res_table] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, excluded, processing_always, reviewFlag, pause_cases)
% [ROITable, sf, rf, sf_i, rf_i, res_table] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, excluded, processing_always, reviewFlag, pause_cases)
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
    pause_cases = 0;
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

sSD = [];
rSD = [];
sSD_i = [];
rSD_i = [];

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

stress_aif_duration = zeros(num,1);
HeartRate_stress = zeros(num, 1);

for n=1:num
       
    ii = n;
    
    stressCase = PerfTable{n+1, stress_column};
    restCase = PerfTable{n+1, rest_column};
    if(ischemia_column<=num_column)
        Category = PerfTable{n+1, ischemia_column};
    else
        Category = 1;
    end
    
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]); 
    disp(['==================================================================']);  
       
    [configName, scannerID, pID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(stressCase);
    [configName, scannerID, pID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(restCase);

    figDir = fullfile(resDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_Figure']) 
    
    if(excluded)
        roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_ROI_excluded'])
    else
        roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_ROI'])
    end
    
    scanInd = [scanInd; n];
    patientID = [patientID; {pID}];
    scanDate = [scanDate; study_dates];
    scanTime = [scanTime; study_time_stress];
    stressHB = [stressHB; PerfTable{n+1, 9}];
    restHB = [restHB; PerfTable{n+1, 13}];
    age = [age; PerfTable{n+1, 14}];
    gender = [gender; PerfTable{n+1, 15}];
    hematocrit = [hematocrit; PerfTable{n+1, hct_column}];

    HCT = PerfTable{n+1, hct_column};
    if(isnan(HCT))
        HCT = 0;
    else
        if(HCT>1)
            HCT = HCT/100;
        end
    end
    
    %% aif related
    try
        [HeartRate_rest(n), aif_cin_Gd_rest, aif_cin_Gd_rest_without_R2Star, aif_cin_all_echo0_signal_rest, aif_cin_all_echo0_signal_after_R2StarCorrection_rest, footTime_rest, peakTime_rest, valleyTime_rest, R2Star_rest, SampledInterval_rest, KiMap, flowMap, EMap, PSMap, VisfMap, VpMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, restCase);
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
    
    [HeartRate_stress(n), aif_cin_Gd_stress, aif_cin_Gd_stress_without_R2Star, aif_cin_all_echo0_signal_stress, aif_cin_all_echo0_signal_after_R2StarCorrection_stress, footTime_stress, peakTime_stress, valleyTime_stress, R2Star_stress, SampledInterval_stress, KiMap, flowMap, EMap, PSMap, VisfMap, VpMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, stressCase);

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
        
    has_roi_1 = isFileExist(fullfile(roiDir, s1_roi)) | isFileExist(fullfile(roiDir, r1_roi));
    has_roi_2 = isFileExist(fullfile(roiDir, s2_roi)) | isFileExist(fullfile(roiDir, r2_roi));
    has_roi_3 = isFileExist(fullfile(roiDir, s3_roi)) | isFileExist(fullfile(roiDir, r3_roi));
    
    if(has_roi_1 | has_roi_2 | has_roi_3)
       disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]);       
        
       two_ROI = 0;

        PerfResult_file = fullfile(roiDir, 'PerfResult.mat');
        if(~processing_always & isFileExist(PerfResult_file))
            tt = load(PerfResult_file);
            v = tt.PerfResult;
        else        
            cd(roiDir)
            v = PerfTable(n+1, :);

            nV = numel(v);

%             s1 = load(fullfile(roiDir, s1_roi));
%             s2 = load(fullfile(roiDir, s2_roi));
%             s3 = load(fullfile(roiDir, s3_roi));
% 
%             r1 = load(fullfile(roiDir, r1_roi));
%             r2 = load(fullfile(roiDir, r2_roi));
%             r3 = load(fullfile(roiDir, r3_roi));

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
            
%             if(~isFileExist(fullfile(roiDir, 'rest.mat')))
%                 copyfile(fullfile(figDir, 'rest.mat'), roiDir);
%             end
%             
%             if(~isFileExist(fullfile(roiDir, 'stress.mat')))
%                 copyfile(fullfile(figDir, 'stress.mat'), roiDir);
%             end
            
            if(has_rest)
                rest = load(fullfile(figDir, 'rest.mat'));
                res_rest = PerformGadgetronRecon_SavedIsmrmrd_ROIValues_OneCase(r1, r2, r3, rest);
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
            res_stress = PerformGadgetronRecon_SavedIsmrmrd_ROIValues_OneCase(s1, s2, s3, stress);
            
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
            
%             v{nV+1} = res_stress.flow;
%             disp(['stress flow : ' num2str(v{nV+1})]);
% 
%             v{nV+2} = res_rest.flow;
%             disp(['rest flow : ' num2str(v{nV+2})]);
% 
%             if(numel(s1.ROI_info_table)==2)
%                 
%                 two_ROI = 1;
% 
%                 v{nV+3} = res_stress.flow_i;
%                 disp(['stress flow, ischemia : ' num2str(v{nV+3})]);
%                 
%                 v{nV+4} = res_rest.flow_i;
%             else
%                 v{nV+3} = [-1 -1 -1];
%                 v{nV+4} = [-1 -1 -1];
%             end

            v{nV+5} = max(stress.aif_stress_baseline_corrected);
            v{nV+6} = max(rest.aif_rest_baseline_corrected);
            v{nV+7} = HCT;
            
            
            if(reviewFlag)
                if(~isempty(s1)) figure; imagescn(stress.flow_stress(:,:,1,end), [0 8], [], [], [], fullfile(roiDir, s1_roi)); PerfColorMap; end
                if(~isempty(s2)) figure; imagescn(stress.flow_stress(:,:,2,end), [0 8], [], [], [], fullfile(roiDir, s2_roi)); PerfColorMap; end
                if(~isempty(s3)) figure; imagescn(stress.flow_stress(:,:,3,end), [0 8], [], [], [], fullfile(roiDir, s3_roi)); PerfColorMap; end

                if(has_rest)
                    if(~isempty(r1)) figure; imagescn(rest.flow_rest(:,:,1,end), [0 8], [], [], [], fullfile(roiDir, r1_roi)); PerfColorMap; end
                    if(~isempty(r2)) figure; imagescn(rest.flow_rest(:,:,2,end), [0 8], [], [], [], fullfile(roiDir, r2_roi)); PerfColorMap; end
                    if(~isempty(r3)) figure; imagescn(rest.flow_rest(:,:,3,end), [0 8], [], [], [], fullfile(roiDir, r3_roi)); PerfColorMap; end
                end
                
                if(pause_cases) 
                    user_in = input('accept cases y or n :');
                    if(user_in=='n')
                        onlyReview = 1;
                        [h_flow_stress, h_flow_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir,  stressCase, restCase, [0 6], onlyReview, resDir);      
                        pause;
                    end
                end
                closeall
            end

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

            sDelay = [sDelay; res_stress.delay];
            rDelay = [rDelay; res_rest.delay];

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
        end
        
%             f1 = roi_statistics(stress.flow_stress(:,:,1), s1.ROI_info_table(1,1));
%             f2 = roi_statistics(stress.flow_stress(:,:,2), s2.ROI_info_table(1,1));
%             f3 = roi_statistics(stress.flow_stress(:,:,3), s3.ROI_info_table(1,1));
% 
%             v{nV+1} = [f1.m f2.m f3.m];
%             disp(['stress flow : ' num2str(v{nV+1})]);
% 
%             f1 = roi_statistics(rest.flow_rest(:,:,1), r1.ROI_info_table(1,1));
%             f2 = roi_statistics(rest.flow_rest(:,:,2), r2.ROI_info_table(1,1));
%             f3 = roi_statistics(rest.flow_rest(:,:,3), r3.ROI_info_table(1,1));
% 
%             v{nV+2} = [f1.m f2.m f3.m];
%             disp(['rest flow : ' num2str(v{nV+2})]);
% 
%             if(numel(s1.ROI_info_table)==2)
%                 
%                 two_ROI = 1;
%                 
%                 if(numel(s1.ROI_info_table)==2)
%                     f1 = roi_statistics(stress.flow_stress(:,:,1), s1.ROI_info_table(2));
%                 else
%                     f1.m = -1;
%                 end
% 
%                 if(numel(s2.ROI_info_table)==2)
%                     f2 = roi_statistics(stress.flow_stress(:,:,2), s2.ROI_info_table(2));
%                 else
%                     f2.m = -1;
%                 end
% 
%                 if(numel(s3.ROI_info_table)==2)
%                     f3 = roi_statistics(stress.flow_stress(:,:,3), s3.ROI_info_table(2));
%                 else
%                     f3.m = -1;
%                 end
% 
%                 v{nV+3} = [f1.m f2.m f3.m];
%                 disp(['stress flow, ischemia : ' num2str(v{nV+3})]);
%             else
%                 v{nV+3} = [];
%                 v{nV+4} = [];
%             end
% 
%             v{nV+5} = max(stress.aif_stress_baseline_corrected);
%             v{nV+6} = max(rest.aif_rest_baseline_corrected);
%         end
    
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
        pause;
    end
end

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
                stress_aif_T2S_peak, stress_aif_T2S_baseline, stress_aif_duration, ...                
                HeartRate_rest, ... 
                rest_aif_peak, rest_aif_peak_no_T2Star, rest_aif_peak_intensity, rest_aif_peak_intensity_no_T2Star, ... 
                rest_aif_valley, rest_aif_valley_no_T2Star, rest_aif_valley_intensity, rest_aif_valley_intensity_no_T2Star, ... 
                rest_aif_T2S_peak, rest_aif_T2S_baseline, rest_aif_duration, ...                 
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
                sf_mean, rf_mean, ... 
                sE_mean, rE_mean, ... 
                sPS_mean, rPS_mean, ...
                sVisf_mean, rVisf_mean, ... 
                sVp_mean, rVp_mean, ... 
                sKi_MF_mean, rKi_MF_mean, ... 
                sKi_Fermi_mean, rKi_Fermi_mean, ... 
                sKi_TwoCompExp_mean, rKi_TwoCompExp_mean, ...
                sKi_BTEX_mean, rKi_BTEX_mean, ...
                sDelay, rDelay);

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
