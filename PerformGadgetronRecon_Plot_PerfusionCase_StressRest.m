
function [h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir, stressCase, restCase, flow_windowing, onlyReview, checkprocessed, baseDir)
% [h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir, stressCase, restCase, flow_windowing, onlyReview, checkprocessed, baseDir)
% PerformGadgetronRecon_Plot_PerfusionCase_StressRest('I:\ReconResults\BARTS', stressCase, restCase, onlyReview)

    if(nargin < 4)
        flow_windowing = [0 8];
    end

    if(nargin < 5)
        onlyReview = 0;
    end

    if(nargin < 6)
        checkprocessed = 1;
    end
    
    if(nargin < 7)
        baseDir = resDir;
    end
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates_stress, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(stressCase);
    [configName, scannerID, patientID, studyID, measurementID, study_dates_rest, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(restCase);

    if strfind(study_dates_stress,'201');
        restDir = fullfile(resDir, study_dates_rest, restCase)
        stressDir = fullfile(resDir, study_dates_stress, stressCase)
    elseif strfind(study_dates_stress,'202');
        restDir = fullfile(resDir, study_dates_rest, restCase)
        stressDir = fullfile(resDir, study_dates_stress, stressCase)
    else
        restDir = fullfile(resDir, [], restCase)
        stressDir = fullfile(resDir, [], stressCase)        
    end
    
    h_flow_stress = -1;
    h_flow_rest = -1;
    has_stress = -1;
    has_rest = -1;
    try
        [h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest_WithInfo(resDir, restDir, stressDir, scannerID, patientID, studyID, study_dates_stress, study_time_stress, study_dates_rest, study_time_rest, flow_windowing, onlyReview, checkprocessed, baseDir);
    catch
    end
    