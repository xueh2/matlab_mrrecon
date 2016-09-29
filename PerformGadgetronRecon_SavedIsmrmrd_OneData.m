
function timeUsed = PerformGadgetronRecon_SavedIsmrmrd_OneData(dataDir, data_name, gt_host, resDir, styleSheet)
% timeUsed = PerformGadgetronRecon_SavedIsmrmrd_OneData(dataDir, data_name, gt_host, resDir, styleSheet)
% timeUsed = PerformGadgetronRecon_SavedIsmrmrd_OneData('I:\KAROLINSKA', 'DB_LGE_MOCO_AVE_OnTheFly_41672_7432041_7432046_590_20160330-113443', 'localhost', 'I:\ReconResults\KAROLINSKA')

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

[key, user] = sshKeyLookup(gt_host);

setenv('GT_HOST', gt_host); setenv('GT_PORT', GT_PORT);

if(nargin<4)
    resDir = dataDir;
end

if(nargin<5)
    styleSheet = 'IsmrmrdParameterMap_Siemens.xsl';
end

getenv('GT_HOST')
getenv('GT_PORT')

GTHome = getenv('GADGETRON_HOME')
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

styleSheetDefault = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl';
styleSheetPerfusionUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xsl';
if ( nargin >= 5 )
    styleSheetDefault = [ '%GADGETRON_DIR%\install\schema/' styleSheet];
end

xmlUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xml';

isPerf = 0;
if(isempty(strfind(data_name, 'Perfusion'))~=1)
    isPerf = 1;
end
    
% ------------------------------------------------------------
  
dataName = fullfile(dataDir, [data_name '.h5']);

[configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(data_name);

try
    noise_mear_id = findNoiseDependencyMeasurementID_SavedIsmrmrd(dataName);
catch
    noise_mear_id = [];
end

if(~isempty(noise_mear_id))

    [names_noise, numNoise] = findFILE(dataDir, ['*' noise_mear_id '*.h5']);

    for kk=numNoise:numNoise
        h5Name = names_noise{kk};  

        finfo = dir(h5Name);

        if(finfo.bytes<5*1024*1024)
            disp(['File size too small - ' num2str(n) ' - ' name]);
            continue;
        end

        command = ['gadgetron_ismrmrd_client -f ' h5Name ' -c default_measurement_dependencies.xml -a %GT_HOST% -p %GT_PORT% ']
        dos(command, '-echo');
    end
end

dstDir = fullfile(resDir, study_dates, data_name);
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

finfo = dir(dataName);

if(isempty(strfind(data_name, 'Perfusion'))~=1)
    if(finfo.bytes<200*1024*1024)
        disp(['File size too small - ' num2str(n) ' - ' name]);
        cd(dataDir)
        delete(fullfile(dstDir, '*.*'));
        rmdir(fullfile(dstDir), 's');
        return;
    end
else
    if(finfo.bytes<20*1024*1024)
        disp(['File size too small - ' num2str(n) ' - ' name]);
        cd(dataDir)
        delete(fullfile(dstDir, '*.*'));
        rmdir(fullfile(dstDir), 's');
        return;
    end
end

%% run the scan
configNameUsed = fullfile(GTConfigFolder, configName);

lenC = length(configName);
currPath = pwd;

lenUsed = 200 - (1+length(currPath));   
if(lenUsed>lenC)
    lenUsed = lenC;
end

configNameShortened = configName(1:lenUsed);

command = ['gadgetron_ismrmrd_client -f ' dataName ' -C ' configNameUsed ' -a %GT_HOST% -p %GT_PORT% -F %OutputFormat% -G ' configNameShortened ' -o ref_' date_suffix '.h5']
tic; dos(command); timeUsed = toc;

if(isPerf)
    debugFolder = 'E:\gtuser\mrprogs\install\DebugOutput';
    try
        rmdir(debugFolder, 's');
    catch
    end

    try
        mkdir(debugFolder);
    catch
    end
end

if(isPerf)
    if(strcmp(getenv('GT_HOST'), 'localhost')==1)
        movefile(debugFolder, dstDir, 'f');
    else
        [key, user] = sshKeyLookup(getenv('GT_HOST'));
        debug_folder = ['/home/' user '/Debug/DebugOutput']
        CopyGadgetronDebugOutputOnRemote(getenv('GT_HOST'), debug_folder, dstDir, 1)
    end
end
    
%% copy the dicom
PerformGadgetronRecon_SavedIsmrmrd_CopyDicom(resDir, data_name, gt_host);

%% open the results
winopen(dstDir)
