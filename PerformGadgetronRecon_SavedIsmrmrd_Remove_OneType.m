
function PerformGadgetronRecon_SavedIsmrmrd_Remove_OneType(dataDir, scan_type, date_start, date_end, resDir)
% PerformGadgetronRecon_SavedIsmrmrd_Remove_OneType(dataDir, scan_type, date_start, date_end, resDir)
% PerformGadgetronRecon_SavedIsmrmrd_Remove_OneType('I:\KAROLINSKA', {'LGE'}, '2016-04-22', '2017-07-11', 'I:\ReconResults\KAROLINSKA')
% PerformGadgetronRecon_SavedIsmrmrd_Remove_OneType('I:\ROYALFREE', {'LGE', 'Perfusion'}, '2016-01-01', '2017-01-01', 'I:\ReconResults\ROYALFREE')
% PerformGadgetronRecon_SavedIsmrmrd_Remove_OneType('I:\BARTS', {'LGE'}, '2016-01-01', '2017-01-01', 'I:\ReconResults\BARTS')

GTHome = getenv('GADGETRON_HOME')
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

% ------------------------------------------------------------

% find data

files = [];

startN = datenum(date_start);
endN = datenum(date_end);

nT = numel(scan_type);

[names, num] = findFILE(dataDir, '*.h5');          
for n=1:num
    
    process = 0;
    for kk=1:nT
        if(~isempty(strfind(names{n}, scan_type{kk})))
            process = 1;
            break;
        end
    end
    
    if(process==0)
        continue;
    end
    
    [pathstr, name, ext] = fileparts(names{n});
    
    % find scanner ID, patient ID, study ID, measurement ID, study date and time
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
       
    if( str2num(measurementID) > 10000 )
        continue;
    end
    tt = datenum(str2num(study_year), str2num(study_month), str2num(study_day));
    
    if (tt<=endN && tt>=startN)
        disp(name);
        files = [files; {name}];
    end
end 

num = numel(files);
for n=1:num

    name = files{n};
    if(~isempty(strfind(files{n}, 'ISMRMRD_Noise_dependency_')))
        continue;
    end
    
    isPerf = 0;
    if(isempty(strfind(name, 'Perfusion'))~=1)
        isPerf = 1;
    end
       
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
        
    dstDir = fullfile(resDir, study_dates, name);
    disp([num2str(n) ' out of ' num2str(num) ' : ' dstDir]);
    try
        if(exist(dstDir) == 7)
            rmdir(dstDir, 's');   
            rmdir([dstDir '_dicom'], 's');   
        end
    catch
        disp(['Failed to remove ' dstDir]);
    end
end
