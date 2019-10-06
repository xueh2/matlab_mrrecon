
function scan_info = PerformGadgetronRecon_Find_ScanStatistics(dataDir, resDir, date_start, date_end)
% scan_info = PerformGadgetronRecon_Find_ScanStatistics(dataDir, resDir, date_start, date_end)
% scan_info = PerformGadgetronRecon_Find_ScanStatistics('\\hl-share\RawMRI\Lab-Kellman\RawData\BARTS', '\\hl-share\RawMRI\Lab-Kellman\ReconResults\BARTS')
% scan_info : [patient_id study_date No. of {LGE, DBLGE, perfusion, binning cine, realtime cine, T2w, T2*}]

if(nargin<3)
    date_start = '2016-01-01';
end

if(nargin<4)
    date_end = '2029-01-01';
end

% ------------------------------------------------------------
  
perf_cases = [];
files = [];

configNames = [];
startN = datenum(date_start);
endN = datenum(date_end);

[subdirs, numdirs] = FindSubDirs(dataDir);
for d=1:numdirs
    
    study_dir = subdirs{d};
    disp(study_dir);
    [names, num] = findFILE(fullfile(dataDir, subdirs{d}), '*.h5');
    
    try
        tt = datenum(str2num(study_dir(1:4)), str2num(study_dir(5:6)), str2num(study_dir(7:8)));
    catch
        continue;
    end
    
    if(~isnumeric(tt) | isempty(tt))
        continue;
    end
    
    if (tt>endN || tt<startN)
        continue;
    end
    
    for n=1:num
%         process = 0;
%         for kk=1:nT
%             if(~isempty(strfind(names{n}, scan_type{kk})))
%                 process = 1;
%                 break;
%             end
%         end
% 
%         if(process==0)
%             continue;
%         end

        [pathstr, name, ext] = fileparts(names{n});

        % find scanner ID, patient ID, study ID, measurement ID, study date and time
        try
            [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
        catch
            continue;
        end
        
        if( str2num(measurementID) > 10000 )
            continue;
        end
        tt = datenum(str2num(study_year), str2num(study_month), str2num(study_day));

        if (tt<=endN && tt>=startN)
            
            finfo = dir(names{n});
            
            try
                if(finfo.bytes>16*1024*1024)
                    % disp(name);
                    files = [files; {name}];
                    configNames = [configNames; {configName}];
                end
            catch
            end
        end
    end
end 

num = numel(files);
disp(['Found in total ' num2str(num) ' stored files ... ']);
disp('======================================================================');
tUsed = [];
ignored = [];
files_processed = [];
study_dates = [];
week_num = [];
sha1_processed = [];
patientID_all = [];
scannerID_all = [];

LGE = [];
DBLGE = [];
Perf = [];
Binning = [];
RTCine = [];
T2W = [];
T2S = [];
RetroCine = [];

num_small_file = 1;

maxW = weeknum(datenum('20161231', 'yyyymmdd'));

for n=1:num

    name = files{n};
    if(~isempty(strfind(files{n}, 'ISMRMRD_Noise_dependency_')))
        continue;
    end
    
    isPerf = 0;
    if(isempty(strfind(name, 'Perfusion'))~=1)
        isPerf = 1;
    end

    [isLGE, isDBLGE, isPerf, isBinning, isRTCine, isT2W, isT2S, isRetroCine] = find_data_type(name);
    [configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);

    dataName = fullfile(dataDir, study_date, [name '.h5']);
               
%     tic
    finfo = dir(dataName);
%     toc

    isPerf = 0;
    if(isempty(strfind(name, 'Perfusion'))~=1)
        isPerf = 1;
        if(finfo.bytes<200*1024*1024)
            disp([num2str(num_small_file) ' - file size too small - ' num2str(n) ' - ' name ' - ' num2str(finfo.bytes/1024/1024) 'Mb']);
            num_small_file = num_small_file + 1;
            continue;
        end
    end
    
    sha1 = 0;
    try
        % find ID
%         dset = ismrmrd.Dataset(dataName);
%         header = ismrmrd.xml.deserialize(dset.readxml());
% 
%         measurementID = header.measurementInformation.measurementID;
%         if(isfield(header, 'subjectInformation'))
%             sha1 = header.subjectInformation.patientID;
%             sha1
%         else
%             sha1 = 0;
%         end
% 
%         dset.close();
    catch
        continue;
    end

    files_processed = [files_processed; {name}];
    
    found_pt = 0;
    for kk=1:numel(patientID_all)
        if( ~isempty(strfind(patientID_all{kk}, [patientID '_' studyID])) )
            found_pt = kk;
            break;
        end
    end
    
    if(found_pt>0)
        if(isLGE) LGE(found_pt) = LGE(found_pt) + 1; end
        if(isDBLGE) DBLGE(found_pt) = DBLGE(found_pt) + 1; end
        if(isPerf) Perf(found_pt) = Perf(found_pt) + 1; end
        if(isBinning) Binning(found_pt) = Binning(found_pt) + 1; end
        if(isRTCine) RTCine(found_pt) = RTCine(found_pt) + 1; end
        if(isT2W) T2W(found_pt) = T2W(found_pt) + 1; end
        if(isT2S) T2S(found_pt) = T2S(found_pt) + 1; end
        if(isRetroCine) RetroCine(found_pt) = RetroCine(found_pt) + 1; end
    else
        patientID_all = [patientID_all; {[patientID '_' studyID]}];
        scannerID_all = [scannerID_all; {scannerID}];
        study_dates = [study_dates; str2double(study_date)];
        sha1_processed = [sha1_processed; {sha1}];
        
        wn = weeknum(datenum(study_date, 'yyyymmdd'));        
        if( wn<maxW & ~isempty(strfind(study_date, '2017')) )
            wn = wn + maxW;
        end
        
        if( wn<maxW & ~isempty(strfind(study_date, '2018')) )
            wn = wn + 2*maxW;
        end

        if( wn<maxW & ~isempty(strfind(study_date, '2019')) )
            wn = wn + 3*maxW;
        end
        
        week_num = [week_num; wn];
        
        if(isLGE) LGE = [LGE; 1]; else LGE = [LGE; 0]; end
        if(isDBLGE) DBLGE = [DBLGE; 1]; else DBLGE = [DBLGE; 0]; end
        if(isPerf) Perf = [Perf; 1]; else Perf = [Perf; 0]; end
        if(isBinning) Binning = [Binning; 1]; else Binning = [Binning; 0]; end
        if(isRTCine) RTCine = [RTCine; 1]; else RTCine = [RTCine; 0]; end
        if(isT2W) T2W = [T2W; 1]; else T2W = [T2W; 0]; end
        if(isT2S) T2S = [T2S; 1]; else T2S = [T2S; 0]; end
         if(isRetroCine) RetroCine = [RetroCine; 1]; else RetroCine = [RetroCine; 0]; end
   end
end

scan_info = table(patientID_all, sha1_processed, study_dates, week_num, scannerID_all, LGE, DBLGE, Perf, Binning, RTCine, T2W, T2S, RetroCine);

end

% {LGE, DBLGE, perfusion, binning cine, realtime cine, T2w, T2*}
function [isLGE, isDBLGE, isPerf, isBinning, isRTCine, isT2W, isT2S, isRetroCine] = find_data_type(name)

    isLGE = 0;
    isDBLGE = 0;
    isPerf = 0;
    isBinning = 0;
    isRTCine = 0;
    isT2W = 0;
    isT2S = 0;
    isRetroCine = 0;
    
    if(~isempty(strfind(name, 'DB_LGE')))
        isDBLGE = 1;
    elseif(~isempty(strfind(name, 'LGE')))
        isLGE = 1;
    elseif(~isempty(strfind(name, 'Perfusion')))
        isPerf = 1;
    elseif(~isempty(strfind(name, 'Binning')))
        isBinning = 1;
    elseif(~isempty(strfind(name, 'Cine_NL')))
        isRTCine = 1;
    elseif(~isempty(strfind(name, 'T2W')))
        isT2W = 1;
    elseif(~isempty(strfind(name, 'T2Star')))
        isT2S = 1;
    elseif(~isempty(strfind(name, 'Retro_Lin_Cine')))
        isRetroCine = 1;
    end
end
