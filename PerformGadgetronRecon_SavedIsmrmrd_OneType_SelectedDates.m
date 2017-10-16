
function [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_OneType_SelectedDates(dataDir, scan_type, dates, gt_host, resDir, cleanRemote, checkProcessed, sendDicom, startRemoteGT, configNamePreset)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_OneType_SelectedDates(dataDir, scan_type, start_date, end_date, gt_host, resDir, cleanRemote, checkProcessed, sendDicom, startRemoteGT, configNamePreset)
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

if(strcmp(gt_host, 'hongkong'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'bermuda'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'gibraltar'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'grenada'))
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

if(nargin<10)
    configNamePreset = [];
end

getenv('GT_HOST')
getenv('GT_PORT')

GTHome = getenv('GADGETRON_HOME')
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

% styleSheetDefault = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl';
% styleSheetPerfusionUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xsl';
% if ( nargin >= 9 )
%     styleSheetDefault = [ '%GADGETRON_DIR%\install\schema/' styleSheet];
% end

xmlUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xml';

if(cleanRemote)
    [key, user] = sshKeyLookup(gt_host);
    gt_command = ['rm -rf /tmp/gadgetron_data/*'];
    command = ['ssh -i ' key ' ' user '@' gt_host ' "' gt_command '"'];
    command
    dos(command, '-echo');    
end

if(~startRemoteGT)
    StartGadgetronOnRemote(gt_host)
end

% ------------------------------------------------------------

% find data

files = [];

configNames = [];
study_dates = [];
study_times = [];

numdirs = numel(dates);

for d=1:numdirs
    
    [names, num] = findFILE(fullfile(dataDir, dates{d}), '*.h5');
    for n=1:num

        [pathstr, name, ext] = fileparts(names{n});

        processed = 0;
        for kk=1:numel(scan_type)
            if( ~isempty(strfind(name, scan_type{kk})) )
                processed = 1;
                break;
            end
        end

        if(~processed)
            continue;
        end

        % find scanner ID, patient ID, study ID, measurement ID, study date and time
        [configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);

        if(strcmp(gt_host, 'localhost')==1)
            [pathstr, configName, ext] = fileparts(configName);
            configName = [configName '_localhost' ext];
        end

        if(~isempty(configNamePreset))
            configName = configNamePreset;
        end

        if( str2num(measurementID) > 10000 )
            continue;
        end
        tt = datenum(str2num(study_year), str2num(study_month), str2num(study_day));

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
disp(['Total ' num2str(num) ' data are found ... ']);
tUsed = [];
ignored = [];
% files_processed = [];
% for n=1:num
%     disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' files{n}]);
%         
%     [tU, ig] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, files{n}, gt_host, resDir, sendDicom, startRemoteGT, styleSheet);
%     if(~isempty(tU))
%         tUsed = [tUsed; {n, files{n}, tU{2}, tU{3}}];
%     end
% end

% [tU, ig] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, files, gt_host, resDir, checkProcessed, sendDicom, startRemoteGT, styleSheet);
[tUsed, ignored, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, files, gt_host, resDir, checkProcessed, sendDicom, startRemoteGT, configNames, []);

