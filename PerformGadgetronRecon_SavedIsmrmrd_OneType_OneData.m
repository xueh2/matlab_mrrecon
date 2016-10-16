
function [tUsed, ignored, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, filename, gt_host, resDir, checkProcessed, sendDicom, startRemoteGT, configName_preset, noise_dat_processed, styleSheet)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, filename, gt_host, resDir, checkProcessed, sendDicom, startRemoteGT, styleSheet)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData('I:\KAROLINSKA', 'xxxx', 'localhost', 'I:\ReconResults\KAROLINSKA')
% setenv('OutputFormat', 'h5')

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

if(strcmp(gt_host, 'andorra'))
    GT_PORT = '9008';
end

setenv('GT_HOST', gt_host); setenv('GT_PORT', GT_PORT);

if(nargin<4)
    resDir = dataDir;
end

if(nargin<5)
    checkProcessed = 1;
end

if(nargin<6)
    sendDicom = 0;
end

if(nargin<7)
    startRemoteGT = 1;
end

if(nargin<8)
    configName_preset = [];
end

if(nargin<9)
    noise_dat_processed = [];
end

if(nargin<10)
    styleSheet = 'IsmrmrdParameterMap_Siemens.xsl';
end

GTHome = getenv('GADGETRON_HOME');
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

styleSheetDefault = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl';
styleSheetPerfusionUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xsl';
if ( nargin >= 4 )
    styleSheetDefault = [ '%GADGETRON_DIR%\install\schema/' styleSheet];
end

xmlUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xml';

files = {filename};
num = numel(files);
tUsed = [];
ignored = [];
files_processed = [];
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
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
    
    if(strcmp(gt_host, 'localhost')==1)
        [pathstr, configName, ext] = fileparts(configName);
        configName = [configName '_localhost' ext];
    end
    
    if(~isempty(configName_preset))
        configName = configName_preset;
    end
    
    dstDir = fullfile(resDir, study_dates, name);
    
    if(checkProcessed)
        if(exist(dstDir)==7)
            goodStatus = 1;

            if(isPerf)
                try
                    a = readGTPlusExportImageSeries_Squeeze(dstDir, 120);
                catch
                    goodStatus = 0;
                end
            end

            [hdrfiles, numhdr] = findFILE(dstDir, '*.img');
                            
            if(~isPerf)
                dicomFolder = fullfile(resDir, study_dates, [name '_dicom']);
                [dcmfiles, numdcm] = findFILE(dicomFolder, '*.dcm');
                if(numdcm==0 || numdcm~=numhdr)
                    goodStatus = 0;
                end
            else
                dicomFolder = fullfile(resDir, study_dates, [name '_dicom']);
                [dcmfiles, numdcm] = findFILE(dicomFolder, '*.dcm');
                if(numdcm==0 || numdcm>numhdr)
                    goodStatus = 0;
                end
            end
            
            if(goodStatus)
                disp([num2str(n) ' out of ' num2str(num) ' - Already Processed : ' name]);
                continue;
            end
        end
    end
    
    finfo = dir(dataName);
       
    isPerf = 0;
    if(isempty(strfind(name, 'Perfusion'))~=1)
        isPerf = 1;
        if(finfo.bytes<200*1024*1024)
            disp(['File size too small - ' num2str(n) ' - ' name]);
%             cd(dataDir)
            ignored = [ignored; {n, name, finfo.bytes/1024}];
            continue;
        end
    else
        if(finfo.bytes<20*1024*1024)
            disp(['File size too small - ' num2str(n) ' - ' name]);
%             cd(dataDir)
            ignored = [ignored; {n, name, finfo.bytes/1024}];
            continue;
        end
    end
    
    % start gadgetron
    if(startRemoteGT)
        if((strcmp(gt_host, 'localhost')==1))
            cd('D:\gtuser\gt_windows_setup')
            command = ['gadgetron -p %GT_PORT% > D:\Temp\record_' GT_PORT '.txt']
            dos([command ' &'])
        else
            [key, user] = sshKeyLookup(gt_host);
            if (~isempty(user) & startRemoteGT)
                StartGadgetronOnRemote(gt_host, GT_PORT);
            end
        end
    end
    
    noise_mear_id = findNoiseDependencyMeasurementID_SavedIsmrmrd(dataName);
    
    if(~isempty(noise_mear_id))
        
        [names_noise, numNoise] = findFILE(dataDir, ['*' noise_mear_id '*.h5']);
        
        if(numNoise>0)
            for kk=numNoise:numNoise
                h5Name = names_noise{kk};  

                finfo = dir(h5Name);

                if(finfo.bytes<5*1024*1024)
                    disp(['File size too small - ' num2str(n) ' - ' name]);
                    continue;
                end

                noise_processed = 0;
                for pp=1:numel(noise_dat_processed)
                    if( isempty(strfind(noise_dat_processed{pp}, h5name)) ~= 1)
                        disp(['Noise processed - ' h5Name]);
                        noise_processed = 1;
                        break;
                    end
                end
                
                if(noise_processed)
                    break;
                end
                
                command = ['gadgetron_ismrmrd_client -f ' h5Name ' -c default_measurement_dependencies.xml -a %GT_HOST% -p %GT_PORT% ']
                dos(command, '-echo');
                
                noise_processed = [noise_processed; {h5Name}];
            end
        end
    end
        
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
    
    %% run the scan
    
%     disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' name]);
    pause(0.5)
    
    files_processed = [files_processed; {name}];
    
    configNameUsed = fullfile(GTConfigFolder, configName);
   
    lenC = length(configName);
    currPath = pwd;

    lenUsed = 200 - (1+length(currPath));   
    if(lenUsed>lenC)
        lenUsed = lenC;
    end

    configNameShortened = configName(1:lenUsed);

    if(isPerf)
        if(strcmp(getenv('GT_HOST'), 'localhost')==1)
            debugFolder = 'D:\gtuser\mrprogs\install\DebugOutput';
            try
                rmdir(debugFolder, 's');
            catch
            end

            try
                mkdir(debugFolder);
            catch
            end
        end
    end
    
    % run the data 
    if(isPerf)
        command = ['gadgetron_ismrmrd_client -f ' dataName ' -C ' configNameUsed ' -a %GT_HOST% -p %GT_PORT% -F %OutputFormat% -G ' configNameShortened ' -o ref_' date_suffix '.h5']
    else
        command = ['gadgetron_ismrmrd_client -f ' dataName ' -C ' configNameUsed ' -a %GT_HOST% -p %GT_PORT% -F hdr -G ' configNameShortened ' -o ref_' date_suffix '.h5']
    end
    tic; dos(command); timeUsed = toc;
           
    if(isPerf)
        if(strcmp(getenv('GT_HOST'), 'localhost')==1)
            movefile(debugFolder, dstDir, 'f');
        else
            [key, user] = sshKeyLookup(getenv('GT_HOST'));
            debug_folder = ['/home/' user '/Debug/DebugOutput']
            CopyGadgetronDebugOutputOnRemote(getenv('GT_HOST'), debug_folder, dstDir, 1)
        end
    end
    
    [tDicom, remoteFolder] = PerformGadgetronRecon_SavedIsmrmrd_CopyDicom(resDir, name, gt_host);
    
    dstDir = fullfile(resDir, study_dates, name);   
    if(strcmp(gt_host, 'localhost')==1)
        if(startRemoteGT)
            command = ['taskkill /F /FI "IMAGENAME eq gadgetron.*"'];
            dos(command)
        end

        copyfile(['D:\Temp\record_' GT_PORT '.txt'], [dstDir '\record_' gt_host '_' GT_PORT '.txt']);
    else
        [key, user] = sshKeyLookup(gt_host);
        if (~isempty(user) & startRemoteGT)
            StopGadgetronOnRemote(gt_host, GT_PORT);
            CopyGadgetronRecordOnRemote(gt_host, GT_PORT, [dstDir '\record_' gt_host '_' GT_PORT '.txt']);
        end
    end

    if(sendDicom)
        dicomServer = 'barbados';
        dicomPort = 11112;
        command = ['D:\gtuser\gt_windows_setup\dcmtk-3.6.0\install_vc14\bin\storescu ' dicomServer '.nhlbi.nih.gov ' num2str(dicomPort) ' ' remoteFolder ' --scan-directories -aec DCM4CHEE --user admin --password admin'];
        tic; dos(command, '-echo'); toc
    end
    
    tUsed = [tUsed; {name, timeUsed, configName}];
end
