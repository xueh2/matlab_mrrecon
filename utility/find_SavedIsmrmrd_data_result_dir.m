function [dataDir, resultDir, study_date] = find_SavedIsmrmrd_data_result_dir(case_id, dataDirBase, resultDirBase)
% from the case id, find the scanner id, data dir, result dir etc.

if(nargin<2)
    dataDirBase = 'T:\RawData';
end

if(nargin<3)
    resultDirBase = 'T:\ReconResults';
end

[configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(case_id);

if(~isempty(strfind(case_id, 'Cloud')))
    dataDir = 'Cloud';
    resultDir = 'Cloud';
    return;
end

scannerID = str2num(scannerID);

switch scannerID
    case {41837, 42110, 66016}
        dataDir = 'BARTS'
        resultDir = 'BARTS'
    case 42363
        dataDir = 'ROYALFREE'
        resultDir = 'RoyalFree'
    case 42170
        dataDir = 'KCL'
        resultDir = 'KCL'
    case 66097
        dataDir = 'LEEDS'
        resultDir = 'LEEDS'
    case 42034
        dataDir = 'LUND'
        resultDir = 'LUND'
    case {41672, 46184}
        dataDir = 'KAROLINSKA'
        resultDir = 'KAROLINSKA'
    case 141303
        dataDir = 'CHENIESMEWS'
        resultDir = 'CHENIESMEWS'
    case 53531
        dataDir = 'CAMPINAS'
        resultDir = 'CAMPINAS'       
    case 42311
        dataDir = 'GEISINGER'
        resultDir = 'GEISINGER'       
end

dataDir = fullfile(dataDirBase, dataDir);
resultDir = fullfile(resultDirBase, resultDir);