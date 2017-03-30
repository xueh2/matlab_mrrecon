
function [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_PatientTabel_Entries(dataDir, patient_table, dat_lines, dat_columns, gt_host, resDir, cleanRemote, checkProcessed, sendDicom, startRemoteGT)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_PatientTabel_Entries(dataDir, patient_table, dat_lines, dat_columns, gt_host, resDir, cleanRemote, checkProcessed, sendDicom, startRemoteGT)
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

if(strcmp(gt_host, 'andorra'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'bermuda'))
    GT_PORT = '9008';
end

setenv('GT_HOST', gt_host); setenv('GT_PORT', GT_PORT);

if(nargin<5)
    resDir = dataDir;
end

if(nargin<6)
    cleanRemote = 0;
end

if(nargin<7)
    checkProcessed = 1;
end

if(nargin<8)
    sendDicom = 0;
end

if(nargin<9)
    startRemoteGT = 1;
end

if(nargin<11)
    styleSheet = 'IsmrmrdParameterMap_Siemens.xsl';
end

getenv('GT_HOST')
getenv('GT_PORT')

GTHome = getenv('GADGETRON_HOME')
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

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

num = numel(dat_lines);
num_dat = numel(dat_columns);

for n=1:num  
    for kk=1:num_dat
        
        [pathstr, name, ext] = fileparts(patient_table{dat_lines(n), dat_columns(kk)});

        % find scanner ID, patient ID, study ID, measurement ID, study date and time
        [configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);

        if(strcmp(gt_host, 'localhost')==1)
            [pathstr, configName, ext] = fileparts(configName);
            configName = [configName '_localhost' ext];
        end

        if( str2num(measurementID) > 10000 )
            continue;
        end

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

num = numel(files);
tUsed = [];
ignored = [];
[tU, ig] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, files, gt_host, resDir, checkProcessed, sendDicom, startRemoteGT, configNames, []);

