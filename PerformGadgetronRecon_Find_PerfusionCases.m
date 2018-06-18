
function [perf_cases, rest_cases, files_un_processed] = PerformGadgetronRecon_Find_PerfusionCases(dataDir, resDir, date_start, date_end, onlyR3)
% [perf_cases, rest_cases, files_un_processed] = PerformGadgetronRecon_Find_PerfusionCases(dataDir, resDir, date_start, date_end)
% [perf_cases, rest_cases, files_un_processed] = PerformGadgetronRecon_Find_PerfusionCases('I:\BARTS', 'I:\ReconResults\BARTS')
% perf_cases in {stress, rest} order

if(nargin<3)
    date_start = '2016-01-01';
end

if(nargin<4)
    date_end = '2019-01-01';
end

if(nargin<5)
    onlyR3 = 0;
end

% ------------------------------------------------------------
  
perf_cases = [];
files = [];

configNames = [];
startN = datenum(date_start);
endN = datenum(date_end);

scan_type = {'Perfusion'};

nT = numel(scan_type);

[subdirs, numdirs] = FindSubDirs(dataDir);
for d=1:numdirs
    subdirs{d}
    
    currN = datenum(num2str(subdirs{d}), 'yyyymmdd');
    if (currN>endN || currN<startN)
        continue;
    end
    
    [names, num] = findFILE(fullfile(dataDir, subdirs{d}), '*.h5');
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

%         if( str2num(measurementID) > 10000 )
%             continue;
%         end
        tt = datenum(str2num(study_year), str2num(study_month), str2num(study_day));

        if (tt<=endN && tt>=startN)
            disp(name);
            files = [files; {name}];
            configNames = [configNames; {configName}];
        end
    end
end 

num = numel(files);
tUsed = [];
ignored = [];
files_processed = [];
study_dates = [];
sha1_processed = [];
patientID_all = [];
failed_cases = [];

num_small_file = 1;

for n=1:num

    name = files{n};
    if(~isempty(strfind(files{n}, 'ISMRMRD_Noise_dependency_')))
        continue;
    end
    
    isPerf = 0;
    if(isempty(strfind(name, 'Perfusion_AIFR3'))~=1)
        isPerf = 1;
    end
    if(isempty(strfind(name, 'Perfusion_AIF_TwoEchoes'))~=1)
        isPerf = 1;
    end
    
    if(~isPerf)
        continue;
    end
    
    [configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);

    dataName = fullfile(dataDir, study_date, [name '.h5']);
        
    dstDir = fullfile(resDir, study_date, name);
    
%     if(exist(dstDir)==7)
%         files_processed = [files_processed; {name}];
%         study_dates = [study_dates; str2double(study_date)];
%         sha1_processed = [sha1_processed; {sha1}];
%         patientID_processed = [patientID_processed; {patientID}];
%     else
        
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
            sha1 = 0;
            try
                % find ID
                dset = ismrmrd.Dataset(dataName);
                header = ismrmrd.xml.deserialize(dset.readxml());

                measurementID = header.measurementInformation.measurementID;
                if(isfield(header, 'subjectInformation'))
                    sha1 = header.subjectInformation.patientID;
                    sha1
                else
                    sha1 = 0;
                end

                dset.close();
            catch
                disp(['failed to read data set : ' dataName])
                failed_cases = [failed_cases; {dataName}];
                continue;
            end
            
            files_processed = [files_processed; {name}];
            study_dates = [study_dates; str2double(study_date)];
            sha1_processed = [sha1_processed; {sha1}];
            patientID_all = [patientID_all; {patientID}];
        end
%     end
end

% sort the file by scan date
[study_dates, ind] = sort(study_dates);
files_processed = files_processed(ind);
sha1_processed = sha1_processed(ind);
patientID_all = patientID_all(ind);

perf_cases = [];
rest_cases = [];

files_un_processed = [];
patientID_processed = [];

numP = 1;
numR = 1;

% try to match the file name to find stress and rest
while (~isempty(files_processed))
    f1 = files_processed{1};
       
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(f1);
    
    sha1 = find_sha1(patientID_all, sha1_processed, patientID);
    
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
        if(max(aif)<0.2)
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
                rest_cases = [rest_cases; numR, {f1}, patientID, study_dates, sha1]; numR = numR + 1;
                patientID_processed = [patientID_processed; {[patientID '_' studyID]}];
            else
                if(numel(other_fs) >= 1 )
                    for pfs=1:numel(other_fs)
                        [configName_2, scannerID_2, patientID_2, studyID_2, measurementID_2, study_dates_2, study_year_2, study_month_2, study_day_2, study_time_2] = parseSavedISMRMRD(other_fs{pfs});

                        if(str2num(study_time)>str2num(study_time_2))
                            perf_cases = [perf_cases; {numP, other_fs{pfs}, f1, patientID, study_dates, sha1}];
                        else
                            perf_cases = [perf_cases; {numP, f1, other_fs{pfs}, patientID, study_dates, sha1}];
                        end

                        numP = numP + 1;

                        patientID_processed = [patientID_processed; {[patientID '_' studyID]}];
                    end
                    
                    if(numel(other_fs) >= 2 )
                        [configName_2, scannerID_2, patientID_2, studyID_2, measurementID_2, study_dates_2, study_year_2, study_month_2, study_day_2, study_time_2] = parseSavedISMRMRD(other_fs{1});
                        [configName_3, scannerID_3, patientID_3, studyID_3, measurementID_3, study_dates_3, study_year_3, study_month_3, study_day_3, study_time_3] = parseSavedISMRMRD(other_fs{2});

                        if(str2num(study_time_3)>str2num(study_time_2))
                            perf_cases = [perf_cases; {numP, other_fs{1}, other_fs{2}, patientID, study_dates, sha1}];
                        else
                            perf_cases = [perf_cases; {numP, other_fs{2}, other_fs{1}, patientID, study_dates, sha1}];
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

end

function sha1 = find_sha1(patientID_picked, sha1_picked, patientID)
    sha1 = 0;
        
    for tt=1:numel(patientID_picked)
        if(strcmp(patientID_picked{tt}, patientID))
            sha1 = sha1_picked{tt};
        end
    end    
end
