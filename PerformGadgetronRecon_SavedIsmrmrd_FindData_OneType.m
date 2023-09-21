
function [files, files_record] = PerformGadgetronRecon_SavedIsmrmrd_FindData_OneType(dataDir, scan_type, start_date, end_date, protocols)
% files = PerformGadgetronRecon_SavedIsmrmrd_FindData_OneType(dataDir, scan_type, start_date, end_date)
% files = PerformGadgetronRecon_SavedIsmrmrd_FindData_OneType('I:\BARTS', {'LGE'}, '2016-01-01', '2017-01-01')
% setenv('OutputFormat', 'h5')

if(nargin<5)
    protocols = [];
end

files_record = [];
% ------------------------------------------------------------

% find data

files = [];

startN = datenum(start_date);
endN = datenum(end_date);

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
    
    disp(['search ' subdirs{d}]);
    
    [names, num] = findFILE(fullfile(dataDir, subdirs{d}), '*.h5');
    for n=1:num

        [pathstr, name, ext] = fileparts(names{n});

        processed = 0;
        for kk=1:numel(scan_type)
            if( ~isempty(strfind(lower(name), lower(scan_type{kk}))) )
                processed = 1;
                break;
            end
        end

        if(~processed)
            continue;
        end

        % find scanner ID, patient ID, study ID, measurement ID, study date and time
        try
            [configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
        catch
            disp(name);
        end
        if( str2num(measurementID) > 10000 )
            continue;
        end
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

if(~isempty(protocols))

    file_names = [];
    study_dates = [];
    study_times = [];
    patientIDs = [];
    studyIDs = [];
    scannerIDs = [];
    prots = [];
    headers = [];
    measurementIDs = [];
    
    % check protocols
    for ii=1:numel(files)                
        
        name = files{ii};
        try
            [configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
        catch
            disp(['---> ' name]);
            continue;
        end
        
        h5name = fullfile(dataDir, study_date, [name '.h5']);
        
        try
            dset = ismrmrd.Dataset(h5name, 'dataset');        
            hdr = ismrmrd.xml.deserialize(dset.readxml);
            prot = hdr.measurementInformation.protocolName;
            dset.close();
        catch
            continue;
        end
        
        disp(['data ' num2str(ii) ' out of ' num2str(numel(files)) ' --- ' name ' - ' prot]);
        
        found_flag = 0;
        for kk=1:numel(protocols)            
            if(~isempty(strfind(lower(prot), lower(protocols{kk}))))
                found_flag = 1;
                break;
            end            
        end
        
        if(found_flag)
            file_names = [file_names; {name}];
            study_dates = [study_dates; study_date];
            study_times = [study_times; study_time];
            prots = [prots; {prot}];
            headers = [headers; {hdr}];
            patientIDs = [patientIDs; {patientID}];
            studyIDs = [studyIDs; {studyID}];
            scannerIDs = [scannerIDs; {scannerID}];
            measurementIDs = [measurementIDs; str2double(measurementID)];
        end
    end
    
    if(numel(file_names)>1)
        files_record = table(file_names, prots, patientIDs, studyIDs, study_dates, study_times, measurementIDs, scannerIDs, headers);
    end
    
end