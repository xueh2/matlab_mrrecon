
function [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd(dataDir, date, gt_host, processNoise, resDir, styleSheet)
% perform gadgetron reconstruction for the whole study
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd(dataDir, date, gt_host, startRemoteGT, processNoise, resDir, styleSheet)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd('I:\KAROLINSKA', '2016-05-04', 'localhost', 0, resDir)

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

if(nargin<4)
    processNoise = 1;
end

if(nargin<5)
    resDir = dataDir;
end

if(nargin<6)
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
currN = datenum(date);

[names, num] = findFILE(dataDir, '*.h5');          
for n=1:num
    
    [pathstr, name, ext] = fileparts(names{n});
    
    % find scanner ID, patient ID, study ID, measurement ID, study date and time
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
       
    if (currN==datenum(str2num(study_year), str2num(study_month), str2num(study_day)))
        disp(name);
        files = [files; {name}];
        configNames = [configNames; {configName}];
    end
end 

num = numel(files);
if(processNoise)
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
end

tUsed = [];
ignored = [];
for n=1:num

    name = files{n};
    if(~isempty(strfind(files{n}, 'ISMRMRD_Noise_dependency_')))
        continue;
    end
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
       
    dstDir = fullfile(resDir, study_dates, name);
    mkdir(dstDir);
    cd(dstDir)

    delete(fullfile(dstDir, 'res*.h5'));
    delete(fullfile(dstDir, 'out*.h5'));
    delete(fullfile(dstDir, '*.xml'));

    delete(fullfile(dstDir, '*.nii'));
    delete(fullfile(dstDir, 'gadgetron_*.hdr'));
    delete(fullfile(dstDir, 'gadgetron_*.img'));
    delete(fullfile(dstDir, 'Generic*.hdr'));
    delete(fullfile(dstDir, 'Generic*.img'));
    delete(fullfile(dstDir, 'GTPrep*.hdr'));
    delete(fullfile(dstDir, 'GTPrep*.img'));
    delete(fullfile(dstDir, 'GT*.hdr'));
    delete(fullfile(dstDir, 'GT*.img'));
    delete(fullfile(dstDir, '*.attrib'));

    delete(fullfile(dstDir, '*.xml'));
        
    dataName = fullfile(dataDir, [name '.h5']);
    
    finfo = dir(dataName);
       
    if(isempty(strfind(name, 'Perfusion'))~=1)
        if(finfo.bytes<200*1024*1024)
            disp(['File size too small - ' num2str(n) ' - ' name]);
            cd(dataDir)
            delete(fullfile(dstDir, '*.*'));
            rmdir(fullfile(dstDir), 's');
            ignored = [ignored; {n, name, finfo.bytes/1024}];
            continue;
        end
    else
        if(finfo.bytes<20*1024*1024)
            disp(['File size too small - ' num2str(n) ' - ' name]);
            cd(dataDir)
            delete(fullfile(dstDir, '*.*'));
            rmdir(fullfile(dstDir), 's');
            ignored = [ignored; {n, name, finfo.bytes/1024}];
            continue;
        end
    end
    
    %% run the scan
    configNameUsed = fullfile(GTConfigFolder, configNames{n});

    configName = configNames{n};
    
    lenC = length(configName);
    currPath = pwd;

    lenUsed = 200 - (1+length(currPath));   
    if(lenUsed>lenC)
        lenUsed = lenC;
    end

    configNameShortened = configName(1:lenUsed);

    command = ['gadgetron_ismrmrd_client -f ' dataName ' -C ' configNameUsed ' -a %GT_HOST% -p %GT_PORT% -F %OutputFormat% -G ' configNameShortened ' -o ref_' date_suffix '.h5']
    tic; dos(command); timeUsed = toc;
            
    tUsed = [tUsed; {n, name, timeUsed, configName}];
end
