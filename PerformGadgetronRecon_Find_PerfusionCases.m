
function [perf_cases, rest_cases, files_un_processed] = PerformGadgetronRecon_Find_PerfusionCases(dataDir, resDir, date_start, date_end)
% [perf_cases, rest_cases, files_un_processed] = PerformGadgetronRecon_Find_PerfusionCases(dataDir, resDir, date_start, date_end)
% [perf_cases, rest_cases, files_un_processed] = PerformGadgetronRecon_Find_PerfusionCases('I:\BARTS', 'I:\ReconResults\BARTS')
% perf_cases in {stress, rest} order

if(nargin<3)
    date_start = '2016-01-01';
end

if(nargin<4)
    date_end = '2018-01-01';
end

% ------------------------------------------------------------
  
perf_cases = [];
files = [];

configNames = [];
startN = datenum(date_start);
endN = datenum(date_end);

scan_type = {'Perfusion'};

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
    end
end 

num = numel(files);
tUsed = [];
ignored = [];
files_processed = [];
study_dates = [];
num_small_file = 1;

for n=1:num

    name = files{n};
    if(~isempty(strfind(files{n}, 'ISMRMRD_Noise_dependency_')))
        continue;
    end
    
    isPerf = 0;
    if(isempty(strfind(name, 'Perfusion'))~=1)
        isPerf = 1;
    end
    
    dataName = fullfile(dataDir, [name '.h5']);
    
    [configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
    
    dstDir = fullfile(resDir, study_date, name);
    
    if(exist(dstDir)==7)
        files_processed = [files_processed; {name}];
        study_dates = [study_dates; str2double(study_date)];
    else
        
        finfo = dir(dataName);
       
        isPerf = 0;
        if(isempty(strfind(name, 'Perfusion'))~=1)
            isPerf = 1;
            if(finfo.bytes<200*1024*1024)
                disp([num2str(num_small_file) ' - file size too small - ' num2str(n) ' - ' name ' - ' num2str(finfo.bytes/1024/1024) 'Mb']);
                num_small_file = num_small_file + 1;
                continue;
            end
        end
    
        if(isPerf)
            files_processed = [files_processed; {name}];
            study_dates = [study_dates; str2double(study_date)];
        end
    end
end

% sort the file by scan date
[study_dates, ind] = sort(study_dates);
files_processed = files_processed(ind);

perf_cases = [];
rest_cases = [];

files_un_processed = [];
patientID_processed = [];

numP = 1;

% try to match the file name to find stress and rest
while (~isempty(files_processed))
    f1 = files_processed{1};
       
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(f1);
    
    dstDir = fullfile(resDir, study_dates, f1);
    
    elem_id = [1];
    other_fs = [];

    has_aif = 1;
    try
        aif = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_cin.hdr'));
        has_aif = 1;
    catch
        aif = 2.0;
    end
    
    if(has_aif)
        if(max(aif)<1.5)
            files_un_processed = [files_un_processed; {f1}];  
        else    
            for ii=2:numel(files_processed)
                f2 = files_processed{ii};

                if( ~isempty(strfind(f2, scannerID)) && ~isempty(strfind(f2, patientID)) && ~isempty(strfind(f2, studyID)) && ~isempty(strfind(f2, study_dates)) )
                    other_fs = [other_fs; {f2}];
                    elem_id = [elem_id; ii];
                end       
            end

            if(numel(other_fs) == 0 )
                % rest
                rest_cases = [rest_cases; {f1}];
                patientID_processed = [patientID_processed; {[patientID '_' studyID]}];
            else
                if(numel(other_fs) >= 1 )
                    for pfs=1:numel(other_fs)
                        [configName_2, scannerID_2, patientID_2, studyID_2, measurementID_2, study_dates_2, study_year_2, study_month_2, study_day_2, study_time_2] = parseSavedISMRMRD(other_fs{pfs});

                        if(str2num(study_time)>str2num(study_time_2))
                            perf_cases = [perf_cases; {numP, other_fs{pfs}, f1}];
                        else
                            perf_cases = [perf_cases; {numP, f1, other_fs{pfs}}];
                        end

                        numP = numP + 1;

                        patientID_processed = [patientID_processed; {[patientID '_' studyID]}];
                    end
                else
                    files_un_processed = [files_un_processed; {f1}];  
                    for ii=1:numel(other_fs)
                        files_un_processed = [files_un_processed; {other_fs{ii}}];  
                    end
                end
            end
        end
    end
    
    file_temp = [];
    for ii=1:numel(files_processed)
        if( isempty(find(elem_id==ii)) == 1 )
            file_temp = [file_temp; {files_processed{ii}}];
        end
    end
    
    files_processed = file_temp;
    disp(['File left : ' num2str(numel(files_processed))]);
end

disp(['Found ' num2str(size(perf_cases, 1)) ' stress/rest cases from ' num2str(numel(patientID_processed)) ' patients ... ']);

