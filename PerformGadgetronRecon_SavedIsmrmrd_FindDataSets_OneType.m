
function [files] = PerformGadgetronRecon_SavedIsmrmrd_FindDataSets_OneType(dataDir, scan_type, start_date, end_date)
% files = PerformGadgetronRecon_SavedIsmrmrd_FindDataSets_OneType(dataDir, scan_type, start_date, end_date)
% files = PerformGadgetronRecon_SavedIsmrmrd_FindDataSets_OneType(dataDir, {'Binning'}, '2018-07-31', '2018-07-31');

if(nargin<3)
    start_date = '2016-01-01';
end

if(nargin<4)
    end_date = '2020-01-01';
end

GTHome = getenv('GADGETRON_HOME')
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

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

