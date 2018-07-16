
function files = PerformGadgetronRecon_SavedIsmrmrd_FindData_OneType(dataDir, scan_type, start_date, end_date)
% files = PerformGadgetronRecon_SavedIsmrmrd_FindData_OneType(dataDir, scan_type, start_date, end_date)
% files = PerformGadgetronRecon_SavedIsmrmrd_FindData_OneType('I:\BARTS', {'LGE'}, '2016-01-01', '2017-01-01')
% setenv('OutputFormat', 'h5')

% ------------------------------------------------------------

% find data

files = [];

startN = datenum(start_date);
endN = datenum(end_date);

study_dates = [];
study_times = [];

[subdirs, numdirs] = FindSubDirs(dataDir);
for d=1:numdirs
    
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

%         if( str2num(measurementID) > 10000 )
%             continue;
%         end
        tt = datenum(str2num(study_year), str2num(study_month), str2num(study_day));

        if (tt<=endN && tt>=startN)
            files = [files; {name}];
            study_dates = [study_dates; str2double(study_date)];
            study_times = [study_times; str2double(study_time)];
        end
    end 
end

% sort the file by scan date
[study_dates, ind] = sort(study_dates);
files = files(ind);
