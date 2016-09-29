
function [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_RunPatientStudy(dataDir, patient_studies, gt_host, resDir, cleanRemote, checkProcessed, sendDicom, startRemoteGT, configName_preset, styleSheet)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_RunPatientStudy(dataDir, patient_studies, gt_host, resDir, cleanRemote, checkProcessed, sendDicom, startRemoteGT, styleSheet)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_RunPatientStudy('I:\KAROLINSKA', patient_studies, 'localhost', 'I:\ReconResults\KAROLINSKA')
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_RunPatientStudy('I:\ROYALFREE',patient_studies, 'samoa', 'I:\ReconResults\ROYALFREE')
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_RunPatientStudy('I:\BARTS', patient_studies, 'samoa', 'I:\ReconResults\BARTS')
% setenv('OutputFormat', 'h5')

if(strcmp(gt_host, 'palau'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'localhost'))
    GT_PORT = '9002';
end

if(strcmp(gt_host, 'denmark'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'samoa'))
    GT_PORT = '9016';
end

if(strcmp(gt_host, 'barbados'))
    GT_PORT = '9008';
end

setenv('GT_HOST', gt_host); setenv('GT_PORT', GT_PORT);

if(nargin<4)
    resDir = dataDir;
end

if(nargin<5)
    cleanRemote = 0;
end

if(nargin<6)
    checkProcessed = 1;
end

if(nargin<7)
    sendDicom = 0;
end

if(nargin<8)
    startRemoteGT = 1;
end

if(nargin<9)
    configName_preset = [];
end

if(nargin<10)
    styleSheet = 'IsmrmrdParameterMap_Siemens.xsl';
end

getenv('GT_HOST')
getenv('GT_PORT')

GTHome = getenv('GADGETRON_HOME')
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

styleSheetDefault = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl';
styleSheetPerfusionUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xsl';
if ( nargin >= 4 )
    styleSheetDefault = [ '%GADGETRON_DIR%\install\schema/' styleSheet];
end

xmlUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xml';

if(cleanRemote)
    [key, user] = sshKeyLookup(gt_host);
    gt_command = ['rm -rf /tmp/gadgetron_data/*'];
    command = ['ssh -i ' key ' ' user '@' gt_host ' "' gt_command '"'];
    command
    dos(command, '-echo');    
end

% ------------------------------------------------------------

% find data

files = [];

configNames = [];
study_dates = [];
study_times = [];

[names, num] = findFILE(dataDir, '*.h5');          
for n=1:num
    
    [pathstr, name, ext] = fileparts(names{n});
    
    % find scanner ID, patient ID, study ID, measurement ID, study date and time
    [configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
    
    if(strcmp(gt_host, 'localhost')==1)
        [pathstr, configName, ext] = fileparts(configName);
        configName = [configName '_localhost' ext];
    end
    
    if( str2num(measurementID) > 10000 )
        continue;
    end
    
    patient_study_str = [scannerID '_' patientID '_' studyID];
    studyFound = 0;
    for kk=1:size(patient_studies, 1)
        if(strfind(patient_studies{kk, 1}, patient_study_str)==1)
            studyFound = 1;
            break;
        end
    end

    if(studyFound)
        files = [files; {name}];
        configNames = [configNames; {configName}];
        study_dates = [study_dates; str2double(study_date)];
        study_times = [study_times; str2double(study_time)];
    end
end 

% sort the file by scan date
[study_dates, ind] = sort(study_dates);
files = files(ind);
configNames = configNames(ind);
study_times = study_times(ind);

% run the cases
num = numel(files);
tUsed = [];
ignored = [];
files_processed = [];
noise_dat_processed = [];
for n=1:num
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' files{n}]);   
    
%     if( isempty(strfind(files{n}, 'Perfusion'))==1 )
%         continue;
%     end
    
    [tU, ig, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, files{n}, gt_host, resDir, checkProcessed, sendDicom, startRemoteGT, configName_preset, noise_dat_processed, styleSheet);
    if(~isempty(tU))
        tUsed = [tUsed; {n, files{n}, tU{2}, tU{3}}];
    end
end
