
function PerformGadgetronRecon_SavedIsmrmrd_NoiseScan(dataDir, date_start, gt_host, date_end)
% perform gadgetron reconstruction for the whole study
% PerformGadgetronRecon_SavedIsmrmrd_NoiseScan(dataDir, date, gt_host)
% PerformGadgetronRecon_SavedIsmrmrd_NoiseScan('I:\KAROLINSKA', '2016-05-04', 'localhost')
% PerformGadgetronRecon_SavedIsmrmrd_NoiseScan('I:\KAROLINSKA', '2016-04-01', 'localhost', '2016-05-11')

if(strcmp(gt_host, 'palau'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'localhost'))
    GT_PORT = '9002';
end

if(strcmp(gt_host, 'denmark'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'samoa'))
    GT_PORT = '9016';
end

if(strcmp(gt_host, 'barbados'))
    GT_PORT = '9008';
end

if nargin < 4
    date_end = date_start;
end

setenv('GT_HOST', gt_host); setenv('GT_PORT', GT_PORT);

getenv('GT_HOST')
getenv('GT_PORT')

GTHome = getenv('GADGETRON_HOME')
% ------------------------------------------------------------

% find data

files = [];
currN = datenum(date_start);
maxN = datenum(date_end);

[names, num] = findFILE(dataDir, '*.h5');          
for n=1:num
    
    if(isempty(strfind(names{n}, 'ISMRMRD_Noise_dependency_')))
        continue;
    end
    
    [pathstr, name, ext] = fileparts(names{n});
    
    % find scanner ID, patient ID, study ID, measurement ID, study date and time
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
       
    tt = datenum(str2num(study_year), str2num(study_month), str2num(study_day));
    
    if (tt>=currN && tt<=maxN)
        disp(name);
        files = [files; {name}];
    end
end 

num = numel(files);
for n=1:num

    name = files{n};
    if(isempty(strfind(files{n}, 'ISMRMRD_Noise_dependency_')))
        continue;
    end

    h5Name = fullfile(dataDir, [name '.h5']);  

    finfo = dir(h5Name);

    if(finfo.bytes<5*1024*1024)
        disp(['File size too small - ' num2str(n) ' - ' name]);
        continue;
    end

    command = ['gadgetron_ismrmrd_client -f ' h5Name ' -c default_measurement_dependencies.xml -a %GT_HOST% -p %GT_PORT% ']
    dos(command, '-echo');
end
