
function PerformGadgetronRecon_SavedIsmrmrd_CopyResults(dataDir, cases, resDir, date_start, date_end, scan_type, dstDir)
% PerformGadgetronRecon_SavedIsmrmrd_CopyResults(dataDir, cases, resDir, dstDir)
% PerformGadgetronRecon_SavedIsmrmrd_CopyResults('I:\BARTS', cases, 'I:\ReconResults\BARTS', date_start, date_end, {'Perfusion'}, 'I:\ReconResults\BARTS\BARTS_SplenicTest')

% ------------------------------------------------------------

% find data

files = [];

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
       
    tt = datenum(str2num(study_year), str2num(study_month), str2num(study_day));
    
    if (tt<=endN && tt>=startN)
        disp(name);
        files = [files; {name}];
        configNames = [configNames; {configName}];
    end
end 

num = numel(files);
for n=1:num

    name = files{n};
    if(~isempty(strfind(files{n}, 'ISMRMRD_Noise_dependency_')))
        continue;
    end
    
    dataName = fullfile(dataDir, [name '.h5']);
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);

    data_dstDir = fullfile(resDir, study_dates, name);
    
    if(~isFileExist(data_dstDir))
        continue;
    end
    
    tt = datenum(str2num(study_year), str2num(study_month), str2num(study_day));
    
    study_hour = str2num(study_time(1:2));
    study_min = str2num(study_time(3:4));
    study_sec = str2num(study_time(5:6));
    
    for d=1:size(cases, 1)
        
        case_date = cases{d, 1};
        case_stress_time = cases{d, 2};
        case_rest_time = cases{d, 3};
                
        ndate = datenum(case_date,'dd/mm/yyyy');
        
        if(ndate==tt)          
            C = strsplit(case_stress_time, ':');
            case_hour = str2num(C{1});
            case_min = str2num(C{2});
            case_sec = str2num(C{3});
            if(study_hour==case_hour & abs(study_min-case_min)<=2 )
                disp(['Copying ' data_dstDir]);
                
                mkdir(fullfile(dstDir, study_dates));
                copyfile([data_dstDir '_dicom'], fullfile(dstDir, study_dates, name), 'f');
            end
            
            C = strsplit(case_rest_time, ':');
            case_hour = str2num(C{1});
            case_min = str2num(C{2});
            case_sec = str2num(C{3});
            if(study_hour==case_hour & abs(study_min-case_min)<=2 )
                disp(['Copying ' data_dstDir]);

                mkdir(fullfile(dstDir, study_dates));
                copyfile([data_dstDir '_dicom'], fullfile(dstDir, study_dates, name), 'f');
            end
        end
    end
end
