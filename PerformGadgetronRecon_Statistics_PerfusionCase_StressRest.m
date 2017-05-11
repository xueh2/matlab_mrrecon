
function [PerfTable, res_table]  = PerformGadgetronRecon_Statistics_PerfusionCase_StressRest(dataDir, resDir, perf_cases, recomputed)
% PerfTable = PerformGadgetronRecon_Statistics_PerfusionCase_StressRest(dataDir, resDir, perf_cases)
% PerfTable = PerformGadgetronRecon_Statistics_PerfusionCase_StressRest('I:\BARTS', 'I:\ReconResults\BARTS', perf_cases)

if(nargin < 4)
    recomputed = 0;
end

num = size(perf_cases, 1);

patientID = [];
rest_time = [];
stress_time = [];
study_date = [];

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

ind_case_all = zeros(num,1);
field_strength_all = zeros(num,1);
stressCase_all = cell(num,1);
restCase_all = cell(num,1);
stress_sha1_all = cell(num,1);
rest_sha1_all = cell(num,1);
scannerID_all = cell(num,1);

PerfTable = {'num', 'PatientID', 'study date', 'stress scan time', 'rest scan time', 'field strength', 'scanner id', ... 
    'stress aif peak Gd', 'stress aif peak Gd without T2* correction', 'stress aif duration (seconds)', 'stress heart rate', ... 
    'rest aif peak Gd', 'rest aif peak Gd without T2* correction', 'rest aif duration (seconds)', 'rest heart rate', ... 
    'patient height', 'patient weight', 'BSA ratio','patient age', 'patient gender', 'patient history/risk factors', 'hematocrit', 'CMR report/findings', 'EF', ...
    'stress data', 'rest data', 'stress SHA1', 'rest SHA1', ... 
    'stress aif peak inten', 'stress aif peak inten without T2*', 'stress aif valley Gd', 'stress aif valley Gd without T2*', 'stress aif valley inten', 'stress aif valley inten without T2*', ... 
    'rest aif peak inten', 'rest aif peak inten without T2*', 'rest aif valley Gd', 'rest aif valley Gd without T2*', 'rest aif valley inten', 'rest aif valley inten without T2*', 'rest T2Star peak', 'rest T2Star baseline', 'stress T2Star peak', 'stress T2Star baseline'};

PerfTable_Wrong = {'num', 'PatientID', 'study date', 'stress scan time', 'rest scan time', 'stress aif peak Gd', 'stress aif duration (seconds)', 'rest aif peak Gd', 'rest aif duration (seconds)', 'stress data', 'rest data'};

ind_case = 1;
for ii=1:num
    
    restCase = perf_cases{ii, 3}
    stressCase = perf_cases{ii, 2}
    
    disp(['Processing ' num2str(ii) ' out of ' num2str(num) ' - ' stressCase ' - ' restCase]);
    
    [configName, scannerID_rest, patientID_rest, studyID_rest, measurementID_rest, study_dates, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(restCase);
    [configName, scannerID, patientID_stress, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(stressCase);
    
    restDir = fullfile(resDir, study_dates, restCase);
    stressDir = fullfile(resDir, study_dates, stressCase);
    
    record_dat_file = fullfile(stressDir, 'Rest_Stress_Record.mat');
    
    patientID = [patientID; {patientID_stress}];
    rest_time = [rest_time; {study_time_rest}];
    stress_time = [stress_time; {study_time_stress}];
    study_date = [study_date; {study_dates}];

    scannerID_all{ii} = scannerID;
    
    if(~recomputed & isFileExist(record_dat_file))
        ls = load(record_dat_file);
        Rest_Stress_Record = ls.Rest_Stress_Record;

        ind_case_all(ii) = Rest_Stress_Record{1};
        field_strength_all(ii) = Rest_Stress_Record{6};
        scannerID_all{ii} = Rest_Stress_Record{7};
        
        ind_item = 8;
        stress_aif_valley(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        stress_aif_valley_no_T2Star = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        stress_aif_duration = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        HeartRate_stress = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;

        rest_aif_peak(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        rest_aif_peak_no_T2Star = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        rest_aif_duration = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        HeartRate_rest = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        
        ind_item = 25;
        stressCase_all(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;;
        restCase_all(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;;
        stress_sha1_all(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;;
        rest_sha1_all(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;;
        
        stress_aif_peak_intensity(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        stress_aif_peak_intensity_no_T2Star(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        stress_aif_valley(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        stress_aif_valley_no_T2Star(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        stress_aif_valley_intensity(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        stress_aif_valley_intensity_no_T2Star(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        rest_aif_peak_intensity(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        rest_aif_peak_intensity_no_T2Star(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        rest_aif_valley(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        rest_aif_valley_no_T2Star(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        rest_aif_valley_intensity(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        rest_aif_valley_intensity_no_T2Star(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        rest_aif_T2S_peak(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        rest_aif_T2S_baseline(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        stress_aif_T2S_peak(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
        stress_aif_T2S_baseline(ii) = Rest_Stress_Record{ind_item}; ind_item = ind_item + 1;
    else
        try
            [HeartRate_rest(ii), aif_cin_Gd_rest, aif_cin_Gd_rest_without_R2Star, aif_cin_all_echo0_signal_rest, aif_cin_all_echo0_signal_after_R2StarCorrection_rest, footTime_rest, peakTime_rest, valleyTime_rest, R2Star_rest, KiMap, flowMap, EMap, PSMap, VisfMap, VpMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, restCase);

            [HeartRate_stress(ii), aif_cin_Gd_stress, aif_cin_Gd_stress_without_R2Star, aif_cin_all_echo0_signal_stress, aif_cin_all_echo0_signal_after_R2StarCorrection_stress, footTime_stress, peakTime_stress, valleyTime_stress, R2Star_stress,KiMap, flowMap, EMap, PSMap, VisfMap, VpMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, stressCase);

            stressCase_all{ii} = stressCase;
            restCase_all{ii} = stressCase;

            rest_aif_peak(ii) = max(aif_cin_Gd_rest);
            rest_aif_peak_no_T2Star(ii) = max(aif_cin_Gd_rest_without_R2Star);
            rest_aif_peak_intensity(ii) = max(aif_cin_all_echo0_signal_after_R2StarCorrection_rest);
            rest_aif_peak_intensity_no_T2Star(ii) = max(aif_cin_all_echo0_signal_rest);
            
            r1 = ceil(valleyTime_rest);
            r2 = ceil(valleyTime_rest-0.5);
            r3 = ceil(valleyTime_rest+0.5);
            
            rest_aif_valley(ii) = (aif_cin_Gd_rest(r1) + aif_cin_Gd_rest(r2) + aif_cin_Gd_rest(r3))/3;
            rest_aif_valley_no_T2Star(ii) = (aif_cin_Gd_rest_without_R2Star(r1) + aif_cin_Gd_rest_without_R2Star(r2) + aif_cin_Gd_rest_without_R2Star(r3))/3;
            rest_aif_valley_intensity(ii) = (aif_cin_all_echo0_signal_after_R2StarCorrection_rest(r1) + aif_cin_all_echo0_signal_after_R2StarCorrection_rest(r2) + aif_cin_all_echo0_signal_after_R2StarCorrection_rest(r3))/3;
            rest_aif_valley_intensity_no_T2Star(ii) = (aif_cin_all_echo0_signal_rest(r1) + aif_cin_all_echo0_signal_rest(r2) + aif_cin_all_echo0_signal_rest(r3))/3;
            
            rest_aif_duration(ii) = (valleyTime_rest - footTime_rest) * 0.5;

            figure; hold on; plot(aif_cin_all_echo0_signal_rest); plot(aif_cin_all_echo0_signal_after_R2StarCorrection_rest, 'r'); hold off
            figure; hold on; plot(aif_cin_Gd_rest); plot(aif_cin_Gd_rest_without_R2Star, 'r');
            
            rest_T2S = 1.0 ./ (R2Star_rest+eps);
            stress_T2S = 1.0 ./ (R2Star_stress+eps);
            
            rest_aif_T2S_peak(ii) = min( [rest_T2S(ceil(peakTime_rest)) rest_T2S(ceil(peakTime_rest+0.5)) rest_T2S(ceil(peakTime_rest-0.5))] );
            stress_aif_T2S_peak(ii) = min( [stress_T2S(ceil(peakTime_stress)) stress_T2S(ceil(peakTime_stress+0.5)) stress_T2S(ceil(peakTime_stress-0.5))] );
            
            v = rest_T2S(3:floor(footTime_rest));
            ind = find(v<100);
            if (~isempty(ind)) rest_aif_T2S_baseline(ii) = mean(v(ind)); end
            
            v = stress_T2S(3:floor(footTime_stress));
            ind = find(v<100);
            if (~isempty(ind)) stress_aif_T2S_baseline(ii) = mean(v(ind)); end
            
            % -----------------------------------
            
            stress_aif_peak(ii) = max(aif_cin_Gd_stress);
            stress_aif_peak_no_T2Star(ii) = max(aif_cin_Gd_stress_without_R2Star);
            stress_aif_peak_intensity(ii) = max(aif_cin_all_echo0_signal_after_R2StarCorrection_stress);
            stress_aif_peak_intensity_no_T2Star(ii) = max(aif_cin_all_echo0_signal_stress);
            
            r1 = ceil(valleyTime_stress);
            r2 = ceil(valleyTime_stress-0.5);
            r3 = ceil(valleyTime_stress+0.5);
            
            stress_aif_valley(ii) = (aif_cin_Gd_stress(r1) + aif_cin_Gd_stress(r2) + aif_cin_Gd_stress(r3))/3;
            stress_aif_valley_no_T2Star(ii) = (aif_cin_Gd_stress_without_R2Star(r1) + aif_cin_Gd_stress_without_R2Star(r2) + aif_cin_Gd_stress_without_R2Star(r3))/3;
            stress_aif_valley_intensity(ii) = (aif_cin_all_echo0_signal_after_R2StarCorrection_stress(r1) + aif_cin_all_echo0_signal_after_R2StarCorrection_stress(r2) + aif_cin_all_echo0_signal_after_R2StarCorrection_stress(r3))/3;
            stress_aif_valley_intensity_no_T2Star(ii) = (aif_cin_all_echo0_signal_stress(r1) + aif_cin_all_echo0_signal_stress(r2) + aif_cin_all_echo0_signal_stress(r3))/3;
            
            stress_aif_duration(ii) = (valleyTime_stress - footTime_stress) * 0.5;

        %     PerfTable = {'num', 'PatientID', 'study date', 'stress scan time', 'rest scan time', ... 
        %     'stress aif peak Gd', 'stress aif peak Gd without T2* correction', 'stress aif duration (seconds)', 'stress heart rate', ... 
        %     'rest aif peak Gd', 'rest aif peak Gd without T2* correction', 'rest aif duration (seconds)', 'rest heart rate', ... 
        %     'patient age', 'patient gender', 'CMR report/findings', 'EF', ...
        %     'stress data', 'rest data'};

            if(rest_aif_peak(ii)>1.2 && stress_aif_peak(ii)>1.2)

                stress_h5 = fullfile(dataDir, study_dates, [stressCase '.h5']);
                rest_h5 = fullfile(dataDir, study_dates, [restCase '.h5']);
                stress_sha1 = [];
                rest_sha1 = [];

                try
                    dset = ismrmrd.Dataset(stress_h5);
                    header = ismrmrd.xml.deserialize(dset.readxml());
                    if(isfield(header, 'subjectInformation'))
                        stress_sha1 = header.subjectInformation.patientID;
                    end
                    dset.close();

                    dset = ismrmrd.Dataset(rest_h5);
                    header = ismrmrd.xml.deserialize(dset.readxml());
                    if(isfield(header, 'subjectInformation'))
                        rest_sha1 = header.subjectInformation.patientID;
                    end
                    
                    field_strength = header.acquisitionSystemInformation.systemFieldStrength_T;
                    
                    dset.close();
                catch
                end

                ind_case_all(ii) = ind_case;
                stress_sha1_all{ii} = stress_sha1;
                rest_sha1_all{ii} = rest_sha1;
                field_strength_all(ii) = field_strength;

                Rest_Stress_Record = {ind_case, patientID_stress, study_dates, study_time_stress, study_time_rest, field_strength, scannerID, ... 
                    stress_aif_peak(ii), stress_aif_peak_no_T2Star(ii), stress_aif_duration(ii), HeartRate_stress(ii), ... 
                    rest_aif_peak(ii), rest_aif_peak_no_T2Star(ii), rest_aif_duration(ii), HeartRate_rest(ii), ...
                    [], [], [], [], [], [], [], [], [], ... 
                    stressCase, restCase, stress_sha1, rest_sha1, ... 
                    stress_aif_peak_intensity(ii), stress_aif_peak_intensity_no_T2Star(ii), ...
                    stress_aif_valley(ii), stress_aif_valley_no_T2Star(ii), stress_aif_valley_intensity(ii), stress_aif_valley_intensity_no_T2Star(ii), ... 
                    rest_aif_peak_intensity(ii), rest_aif_peak_intensity_no_T2Star(ii), ...
                    rest_aif_valley(ii), rest_aif_valley_no_T2Star(ii), rest_aif_valley_intensity(ii), rest_aif_valley_intensity_no_T2Star(ii), ... 
                    rest_aif_T2S_peak(ii), rest_aif_T2S_baseline(ii), stress_aif_T2S_peak(ii), stress_aif_T2S_baseline(ii)};

                save(fullfile(stressDir, 'Rest_Stress_Record.mat'), 'Rest_Stress_Record');
            else
                disp(['AIF is too low - ' stressCase , ' - ' restCase]);
                
                % pause;
                PerfTable_Wrong = [PerfTable_Wrong; {ii, patientID, study_dates, study_time_stress, study_time_rest, stress_aif_peak(ii), stress_aif_duration(ii), rest_aif_peak(ii), rest_aif_duration(ii), stressCase, restCase}];
            end

            PerfTable = [PerfTable; Rest_Stress_Record];

            ind_case = ind_case + 1;
        catch
            disp('!!! Error happend ... !!!');
        end
    end
    
    closeall
end

PerfTable_Wrong

res_table = table(ind_case_all, patientID, study_date, stress_time, rest_time, field_strength_all, scannerID_all, ... 
    stress_aif_peak, stress_aif_duration, stress_aif_peak_no_T2Star, HeartRate_stress, ...
    rest_aif_peak, rest_aif_duration, rest_aif_peak_no_T2Star, HeartRate_rest, ...
    stressCase_all, restCase_all, stress_sha1_all, rest_sha1_all, ... 
    stress_aif_peak_intensity, stress_aif_peak_intensity_no_T2Star, ...
    stress_aif_valley, stress_aif_valley_no_T2Star, stress_aif_valley_intensity, stress_aif_valley_intensity_no_T2Star, ... 
    rest_aif_peak_intensity, rest_aif_peak_intensity_no_T2Star, ...
    rest_aif_valley, rest_aif_valley_no_T2Star, rest_aif_valley_intensity, rest_aif_valley_intensity_no_T2Star, ... 
    rest_aif_T2S_peak, rest_aif_T2S_baseline, stress_aif_T2S_peak, stress_aif_T2S_baseline);

valid_ind = find(rest_aif_peak(:)>1.5 & stress_aif_peak(:)>1.5);

disp(['Aif peak Gd (mmol/L),              rest   - ' num2str(mean(rest_aif_peak(valid_ind))) ' +/- ' num2str(std(rest_aif_peak(valid_ind)))]);
disp(['Aif peak Gd no T2* (mmol/L),       rest   - ' num2str(mean(rest_aif_peak_no_T2Star(valid_ind))) ' +/- ' num2str(std(rest_aif_peak_no_T2Star(valid_ind)))]);
disp(['Aif duration (seconds),            rest   - ' num2str(mean(rest_aif_duration(valid_ind))) ' +/- ' num2str(std(rest_aif_duration(valid_ind)))]);
disp(['Heart rate,                        rest   - ' num2str(mean(HeartRate_rest(valid_ind))) ' +/- ' num2str(std(HeartRate_rest(valid_ind)))]);

disp(['Aif peak Gd (mmol/L),              stress - ' num2str(mean(stress_aif_peak(valid_ind))) ' +/- ' num2str(std(stress_aif_peak(valid_ind)))]);
disp(['Aif peak Gd no T2* (mmol/L),       stress - ' num2str(mean(stress_aif_peak_no_T2Star(valid_ind))) ' +/- ' num2str(std(stress_aif_peak_no_T2Star(valid_ind)))]);
disp(['Aif duration (seconds),            stress - ' num2str(mean(stress_aif_duration(valid_ind))) ' +/- ' num2str(std(stress_aif_duration(valid_ind)))]);
disp(['Heart rate,                        stress - ' num2str(mean(HeartRate_stress(valid_ind))) ' +/- ' num2str(std(HeartRate_stress(valid_ind)))]);
