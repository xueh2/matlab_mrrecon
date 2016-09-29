
function patient_studies = PerformGadgetronRecon_SavedIsmrmrd_OneType_FindPatietnStudy(dataDir, scan_type, date_start, date_end)
% patient_studies = PerformGadgetronRecon_SavedIsmrmrd_OneType_FindPatietnStudy(dataDir, scan_type, date_start, date_end)
% patient_studies = PerformGadgetronRecon_SavedIsmrmrd_OneType_FindPatietnStudy('I:\KAROLINSKA', {'LGE'}, '2016-04-22', '2016-05-11')
% patient_studies = PerformGadgetronRecon_SavedIsmrmrd_OneType_FindPatietnStudy('I:\ROYALFREE', {'LGE', 'Perfusion'}, '2016-01-01', '2017-01-01')
% patient_studies = PerformGadgetronRecon_SavedIsmrmrd_OneType_FindPatietnStudy('I:\BARTS', {'Perfusion'}, '2016-01-01', '2017-01-01')
% setenv('OutputFormat', 'h5')

GTHome = getenv('GADGETRON_HOME')
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

% ------------------------------------------------------------

% find data

files = [];
patient_studies = [];

configNames = [];
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
        configNames = [configNames; {configName}];
        
        patient_study_str = [scannerID '_' patientID '_' studyID];
        studyFound = 0;
        for kk=1:numel(patient_studies)
            if(strfind(patient_studies{kk}, patient_study_str)==1)
                studyFound = 1;
                break;
            end
        end
        
        if(~studyFound)
            patient_studies = [patient_studies; {patient_study_str study_dates}];
        end
    end
end 
