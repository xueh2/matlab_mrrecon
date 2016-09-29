
function PerformGadgetronRecon_SavedIsmrmrd_CopyDicomAll(dataDir, scan_type, date_start, date_end, gt_host, resDir)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_CopyDicomAll(dataDir, scan_type, date_start, date_end, gt_host, resDir)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_CopyDicomAll('I:\KAROLINSKA', {'LGE'}, '2016-04-22', '2016-05-11', 'localhost', 'I:\ReconResults\KAROLINSKA')
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_CopyDicomAll('I:\ROYALFREE', {'LGE', 'Perfusion'}, '2016-01-01', '2017-01-01', 'samoa', 'I:\ReconResults\ROYALFREE')
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_CopyDicomAll('I:\BARTS', {'Perfusion'}, '2016-01-01', '2017-01-01', 'samoa', 'I:\ReconResults\BARTS')

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

setenv('GT_HOST', gt_host); setenv('GT_PORT', GT_PORT);

if(nargin<6)
    resDir = dataDir;
end

if(nargin<7)
    styleSheet = 'IsmrmrdParameterMap_Siemens.xsl';
end

getenv('GT_HOST')
getenv('GT_PORT')

GTHome = getenv('GADGETRON_HOME')
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

styleSheetDefault = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl';
styleSheetPerfusionUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xsl';
if ( nargin >= 4 )
    styleSheetDefault = [ '%GADGETRON_DIR%\install\schema/' styleSheet];
end

xmlUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xml';
    
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

for n=1:num

    name = files{n};
    if(~isempty(strfind(files{n}, 'ISMRMRD_Noise_dependency_')))
        continue;
    end
    
    dataName = fullfile(dataDir, [name '.h5']);
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
              
    finfo = dir(dataName);
       
    if(isempty(strfind(name, 'Perfusion'))~=1)
        if(finfo.bytes<200*1024*1024)
            disp(['File size too small - ' num2str(n) ' - ' name]);
            continue;
        end
    else
        if(finfo.bytes<20*1024*1024)
            disp(['File size too small - ' num2str(n) ' - ' name]);
            continue;
        end
    end
       
    PerformGadgetronRecon_SavedIsmrmrd_CopyDicom(resDir, name, gt_host);
end
