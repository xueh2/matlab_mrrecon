
function PerfTable = PerformGadgetronRecon_Statistics_UpdatePatientTable_StressRest(resDir, perf_cases)
% PerfTable = PerformGadgetronRecon_Statistics_UpdatePatientTable_StressRest(resDir, PerfTable)
% PerfTable = PerformGadgetronRecon_Statistics_UpdatePatientTable_StressRest('D:\data\ut\NewData\PaperResults\KAROLINSKA_Area_Aif_recorded', PerfTable)

num = size(perf_cases, 1);

patientID = [];
rest_time = [];
stress_time = [];

rest_aif_peak = zeros(num,1);
rest_aif_peak_no_T2Star = zeros(num,1);
rest_aif_duration = zeros(num,1);
HeartRate_rest = zeros(num, 1);

stress_aif_peak = zeros(num,1);
stress_aif_peak_no_T2Star = zeros(num,1);
stress_aif_duration = zeros(num,1);
HeartRate_stress = zeros(num, 1);

sampleinterval = 0.5;
sigmas = [1.6 4.0 5.3];
sigmaMeasure = 0.1;
thresGradient = 0.5;
plotFlag = 0;

PerfTable = {'num', 'PatientID', 'study date', 'stress scan time', 'rest scan time', ... 
    'stress aif peak Gd', 'stress aif peak Gd without T2* correction', 'stress aif duration (seconds)', 'stress heart rate', ... 
    'rest aif peak Gd', 'rest aif peak Gd without T2* correction', 'rest aif duration (seconds)', 'rest heart rate', ... 
    'patient age', 'patient gender', 'CMR report/findings', 'EF', ...
    'stress data', 'rest data'};

PerfTable_Wrong = {'num', 'PatientID', 'study date', 'stress scan time', 'rest scan time', 'stress aif peak Gd', 'stress aif duration (seconds)', 'rest aif peak Gd', 'rest aif duration (seconds)', 'stress data', 'rest data'};

ind = 1;
for ii=1:num-1
    
    restCase = perf_cases{ii+1, 19}
    stressCase = perf_cases{ii+1, 18}
    
    disp([num2str(ii) ' out of ' num2str(num-1) ' - ' stressCase ' - ' restCase]);
    
    [configName, scannerID_rest, patientID_rest, studyID_rest, measurementID_rest, study_dates, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(restCase);
    [configName, scannerID, patientID_stress, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(stressCase);
    
    restDir = fullfile(resDir, study_dates, restCase);
    stressDir = fullfile(resDir, study_dates, stressCase);
    
    [HeartRate_rest(ii), aif_cin_Gd_rest, aif_cin_Gd_without_R2Star_baseline_corrected_rest, footTime_rest, peakTime_rest, valleyTime_rest, KiMap, flowMap, EMap, PSMap, VisfMap, VpMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, restCase);
    
    [HeartRate_stress(ii), aif_cin_Gd_stress, aif_cin_Gd_without_R2Star_baseline_corrected_stress, footTime_stress, peakTime_stress, valleyTime_stress, KiMap, flowMap, EMap, PSMap, VisfMap, VpMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, stressCase);
       
    patientID = [patientID; {patientID_stress}];
    rest_time = [rest_time; {study_time_rest}];
    stress_time = [stress_time; {study_time_stress}];
       
    rest_aif_peak(ii) = max(aif_cin_Gd_rest);
    rest_aif_peak_no_T2Star(ii) = max(aif_cin_Gd_without_R2Star_baseline_corrected_rest);
    rest_aif_duration(ii) = (valleyTime_rest - footTime_rest) * 0.5;
       
    stress_aif_peak(ii) = max(aif_cin_Gd_stress);
    stress_aif_peak_no_T2Star(ii) = max(aif_cin_Gd_without_R2Star_baseline_corrected_stress);
    stress_aif_duration(ii) = (valleyTime_stress - footTime_stress) * 0.5;
    
%     PerfTable = {'num', 'PatientID', 'study date', 'stress scan time', 'rest scan time', ... 
%     'stress aif peak Gd', 'stress aif peak Gd without T2* correction', 'stress aif duration (seconds)', 'stress heart rate', ... 
%     'rest aif peak Gd', 'rest aif peak Gd without T2* correction', 'rest aif duration (seconds)', 'rest heart rate', ... 
%     'patient age', 'patient gender', 'CMR report/findings', 'EF', ...
%     'stress data', 'rest data'};

    if(rest_aif_peak(ii)>1.5 && stress_aif_peak(ii)>1.5)
        
        PerfTable = [PerfTable; {ind, patientID_stress, study_dates, study_time_stress, study_time_rest, ... 
            stress_aif_peak(ii), stress_aif_peak_no_T2Star(ii), stress_aif_duration(ii), HeartRate_stress(ii), ... 
            rest_aif_peak(ii), rest_aif_peak_no_T2Star(ii), rest_aif_duration(ii), HeartRate_rest(ii), ...
            [], [], [], [], ... 
            stressCase, restCase}];
        
        ind = ind + 1;
    else
        disp(['AIF is too low - ' stressCase , ' - ' restCase]);
        pause;
        PerfTable_Wrong = [PerfTable_Wrong; {ii, patientID, study_dates, study_time_stress, study_time_rest, stress_aif_peak(ii), stress_aif_duration(ii), rest_aif_peak(ii), rest_aif_duration(ii), stressCase, restCase}];
    end
    
    closeall
end

PerfTable_Wrong
    
valid_ind = find(rest_aif_peak(:)>1.5 & stress_aif_peak(:)>1.5);

disp(['Aif peak Gd (mmol/L),              rest   - ' num2str(mean(rest_aif_peak(valid_ind))) ' +/- ' num2str(std(rest_aif_peak(valid_ind)))]);
disp(['Aif peak Gd no T2* (mmol/L),       rest   - ' num2str(mean(rest_aif_peak_no_T2Star(valid_ind))) ' +/- ' num2str(std(rest_aif_peak_no_T2Star(valid_ind)))]);
disp(['Aif duration (seconds),            rest   - ' num2str(mean(rest_aif_duration(valid_ind))) ' +/- ' num2str(std(rest_aif_duration(valid_ind)))]);
disp(['Heart rate,                        rest   - ' num2str(mean(HeartRate_rest(valid_ind))) ' +/- ' num2str(std(HeartRate_rest(valid_ind)))]);

disp(['Aif peak Gd (mmol/L),              stress - ' num2str(mean(stress_aif_peak(valid_ind))) ' +/- ' num2str(std(stress_aif_peak(valid_ind)))]);
disp(['Aif peak Gd no T2* (mmol/L),       stress - ' num2str(mean(stress_aif_peak_no_T2Star(valid_ind))) ' +/- ' num2str(std(stress_aif_peak_no_T2Star(valid_ind)))]);
disp(['Aif duration (seconds),            stress - ' num2str(mean(stress_aif_duration(valid_ind))) ' +/- ' num2str(std(stress_aif_duration(valid_ind)))]);
disp(['Heart rate,                        stress - ' num2str(mean(HeartRate_stress(valid_ind))) ' +/- ' num2str(std(HeartRate_stress(valid_ind)))]);
