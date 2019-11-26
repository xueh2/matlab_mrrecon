
function [files, tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_OneType(dataDir, scan_type, start_date, end_date, gt_host, resDir, ... 
    cleanRemote, checkProcessed, sendDicom, startRemoteGT, copy_debug_output, copy_dicom_output, configNamePreset, gt_port_used)
% [files, tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_OneType(dataDir, scan_type, start_date, end_date, gt_host, resDir, cleanRemote, checkProcessed, sendDicom, startRemoteGT, copy_debug_output, copy_dicom_output, configNamePreset, gt_port_used)
% setenv('OutputFormat', 'h5')

GT_PORT = gtPortLookup(gt_host);

% if(strcmp(gt_host, 'palau'))
%     GT_PORT = '9008';
% end
% 
% if(strcmp(gt_host, 'localhost'))
%     GT_PORT = '9002';
% end
% 
% if(strcmp(gt_host, 'denmark'))
%     GT_PORT = '9008';
% end
% 
% if(strcmp(gt_host, 'samoa'))
%     GT_PORT = '9016';
% end
% 
% if(strcmp(gt_host, 'barbados'))
%     GT_PORT = '9008';
% end
% 
% if(strcmp(gt_host, 'andorra'))
%     GT_PORT = '9008';
% end
% 
% if(strcmp(gt_host, 'hongkong'))
%     GT_PORT = '9008';
% end
% 
% if(strcmp(gt_host, 'bermuda'))
%     GT_PORT = '9008';
% end
% 
% if(strcmp(gt_host, 'gibraltar'))
%     GT_PORT = '9008';
% end
% 
if(strcmp(gt_host, 'gt1'))
    gt_host = '137.187.135.97';
end

setenv('GT_HOST', gt_host); 

if(nargin<6)
    resDir = dataDir;
end

if(nargin<7)
    cleanRemote = 0;
end

if(nargin<8)
    checkProcessed = 1;
end

if(nargin<9)
    sendDicom = 0;
end

if(nargin<10)
    startRemoteGT = 1;
end

if(nargin<11)
    copy_debug_output = 1;
end

if(nargin<12)
    copy_dicom_output = 1;
end

if(nargin<13)
    configNamePreset = [];
end

if(nargin<14)
    gt_port_used = GT_PORT;
end

GT_PORT = gt_port_used;
setenv('GT_PORT', GT_PORT);

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
    command = ['ssh -o ConnectTimeout=10 -o StrictHostKeyChecking=no -o TCPKeepAlive=yes -o ServerAliveInterval=15 -o ServerAliveCountMax=3 ' user '@' gt_host ' "' gt_command '"'];
    command
    dos(command, '-echo');    
end

if(startRemoteGT)
    StartGadgetronOnRemote(gt_host)
end

% ------------------------------------------------------------

% find data

files = [];

if(isempty(strfind(start_date, '-')))
    startN = datenum(start_date, 'yyyymmdd');
    endN = datenum(end_date, 'yyyymmdd');
else
    startN = datenum(start_date);
    endN = datenum(end_date);
end

configNames = [];
study_dates = [];
study_times = [];

[subdirs, numdirs] = FindSubDirs(dataDir);
for d=1:numdirs
    
    try
        currN = datenum(subdirs{d}, 'yyyymmdd');

        if(currN<startN | currN>endN)
            continue;
        end
    catch
        continue;
    end
        
    curr_dir = subdirs{d};
    disp(['Search ' subdirs{d} ' - ' num2str(d) ' out of ' num2str(numdirs)])
    tt = datenum(str2num(curr_dir(1:4)), str2num(curr_dir(5:6)), str2num(curr_dir(7:8)));
    if(tt<startN | tt>endN)
        continue;
    end
    
    [names, num] = findFILE(fullfile(dataDir, subdirs{d}), '*.h5');
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

%         if(strcmp(gt_host, 'localhost')==1)
%             [pathstr, configName, ext] = fileparts(configName);
%             configName = [configName '_localhost' ext];
%         end

        if(~isempty(configNamePreset))
            configName = configNamePreset;
        end

%         if( str2num(measurementID) > 10000 )
%             continue;
%         end
        tt = datenum(str2num(study_year), str2num(study_month), str2num(study_day));

        if (tt<=endN && tt>=startN)
            files = [files; {name}];
            configNames = [configNames; {configName}];
            study_dates = [study_dates; str2double(study_date)];
            study_times = [study_times; str2double(study_time)];
        end
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
[tUsed, ignored, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, files, gt_host, resDir, checkProcessed, sendDicom, startRemoteGT, configNames, [], GT_PORT, copy_debug_output, copy_dicom_output);

