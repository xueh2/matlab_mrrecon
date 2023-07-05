% 
function [tUsed, ignored, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData_v2(dataDir, filename, ...
    gt_host, resDir, resDir_local, checkProcessed, delete_old_res, ...
    configName_preset, noise_dat_processed, gt_port, copy_debug_output, copy_dicom_output, pre_set_debug_folder)

% [tUsed, ignored, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData_v2(dataDir, filename, ...
%     gt_host, resDir, resDir_local, checkProcessed, delete_old_res, ...
%     configName_preset, noise_dat_processed, gt_port, copy_debug_output, copy_dicom_output, pre_set_debug_folder)

if(strcmp(gt_host, 'gt1'))
    gt_host = '137.187.135.169';
end
if(strcmp(gt_host, 'gt2'))
    gt_host = '137.187.135.238';
end
if(strcmp(gt_host, 'gt3'))
    gt_host = '137.187.135.194';
end
if(strcmp(gt_host, 'beast'))
    gt_host = '137.187.135.157';
end

GT_PORT = gtPortLookup(gt_host);

setenv('GT_HOST', gt_host); 

output_format = getenv('OutputFormat');

if(isempty(output_format))
    output_format = 'hdr';
end

if(nargin<4)
    resDir = dataDir;
end

if(nargin<5)
    resDir_local = [];
end

if(nargin<6)
    checkProcessed = 1;
end

if(nargin<7)
    delete_old_res = 1;
end

if(nargin<8)
    configName_preset = [];
end

if(nargin<9)
    noise_dat_processed = [];
end

if(nargin<10)
    gt_port = GT_PORT;
end
if(isempty(gt_port))
    gt_port = GT_PORT;
end

if(nargin<11)
    copy_debug_output = 0;
end

if(nargin<12)
    copy_dicom_output = 1;
end

if(nargin<13)
    pre_set_debug_folder = [];
end

setenv('GT_PORT', gt_port);
GT_PORT = gt_port;

GTHome = getenv('GADGETRON_HOME');
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

xmlUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xml';

[key, user] = sshKeyLookup(gt_host);

files = filename;
num = numel(files);
tUsed = [];
ignored = [];
files_processed = [];
noise_id_processed = [];

is_remote_computer = IsRemoteComputer(gt_host);

[keyfile, user] = sshKeyLookup(gt_host);

for n=1:num

    name = files{n};
    if(isnan(name))
        continue;
    end
    if(~isempty(strfind(files{n}, 'ISMRMRD_Noise_dependency_')))
        continue;
    end
    
    isPerf = 0;
    if(isempty(strfind(name, 'Perfusion'))~=1)
        isPerf = 1;
    end
       
    try
        [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
    catch
        disp(['error in parsing ' name])
        continue;
    end
    
%     if(str2num(measurementID)>300000)
%         continue;
%     end
    
    dataName = fullfile(dataDir, study_dates, [name '.h5']);
   
    if(~isempty(configName_preset{1}))
        if(n>numel(configName_preset))
            configName = configName_preset{end};
        else
            configName = configName_preset{n};
        end
    end

    if(iscell(configName))
        configName = configName{1};
    end
    
    disp(['Processing ' num2str(n) ' out of ' num2str(num) ' cases : ' name ' - ' configName]);
    disp(['========================================================================================']);
    
    dstDir = fullfile(resDir, study_dates, name);
    
    if(checkProcessed)
        tstart = tic;
        if(exist(dstDir)==7)
            goodStatus = 1;

            if(isPerf)
                aif_cin_Gd_baseline_corrected = 2;
                try
                    tic
                    aif_cin = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_cin.hdr'));
                    aif_cin_Gd = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
                    aif_cin_Gd_without_R2Star = analyze75read(fullfile(dstDir, 'DebugOutput', 'cin_all_echo0_without_R2Star_LUTCorrection.hdr'));
                    aif_cin_Gd_baseline_corrected = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_cin_echo0_all_signal_baseline_corrected.hdr'));
                    
                    if(max(aif_cin_Gd_baseline_corrected)<1.2)
                        disp('max(aif_cin_Gd_baseline_corrected)<1.2')
                        goodStatus = 1;
                    else                    
                        aif_cin_all_echo0_signal = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_cin_all_echo0_signal.hdr'));
                        aif_cin_all_echo1_signal = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_cin_all_echo1_signal.hdr'));
                        aif_cin_all_echo0_signal_after_R2StarCorrection = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_cin_all_echo0_signal_after_R2StarCorrection.hdr'));
                        aif_cin_all_echo0_OverPD_after_R2StarCorrection = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_cin_all_echo0_OverPD_after_R2StarCorrection.hdr'));
                        aif_cin_all_R2Star = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_cin_all_R2Star.hdr'));
                        aif_cin_all_R2Star_SLEP = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_cin_all_R2Star_SLEP.hdr'));
                        cin_used_all_for_computeFlowMap = analyze75read(fullfile(dstDir, 'DebugOutput', 'cin_used_all_for_computeFlowMap.hdr'));
                        aif_PD = analyze75read(fullfile(dstDir, 'DebugOutput', 'aifPD_for_TwoEcho_T2StartCorrection_0.hdr'));
                        aif_mask = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_LV_mask_for_TwoEcho_T2StartCorrection_0.hdr'));
                        aif_mask_final = analyze75read(fullfile(dstDir, 'DebugOutput', 'AifLVMask_after_Picking.hdr'));
                        PDMaskIntensity = analyze75read(fullfile(dstDir, 'DebugOutput', 'PDMaskIntensity.hdr'));
                        pdPicked = analyze75read(fullfile(dstDir, 'DebugOutput', 'pdPicked.hdr'));
                        aifPicked = analyze75read(fullfile(dstDir, 'DebugOutput', 'aifPicked.hdr'));
                        aifMaskIntensity = analyze75read(fullfile(dstDir, 'DebugOutput', 'aifMaskIntensity.hdr'));
                        aif_LUT = analyze75read(fullfile(dstDir, 'DebugOutput', 'aif_cin_LUT_Valid.hdr'));
                    end
                    disp(['Load hdr to check : ' num2str(toc)]);
                catch
                    goodStatus = 0;
                end
            end

            tic
            if(~isPerf)
                hdrFolder = fullfile(resDir, study_dates, [name]);
                dicomFolder = fullfile(resDir, study_dates, [name '_dicom']);
                if(~isFileExist(dicomFolder))
                    goodStatus = 0;
                    if(isFileExist(hdrFolder))
                        if(strcmp(output_format, 'hdr'))
                            try
                                [dcmfiles, numdcm] = findFILE(dstDir, '*.img');
                            catch
                                numdcm = 0;
                            end
                        else
                            try
                                [dcmfiles, numdcm] = findFILE(dstDir, 'ref*.h5');
                            catch
                                numdcm = 0;
                            end
                        end
                        
                        if(numdcm==0)
                            goodStatus = 0;
                        else
                            goodStatus = 1;
                        end
                        
                        if(copy_dicom_output)
                            goodStatus = 0;
                        end
                    end
                else
                    if(strcmp(output_format, 'hdr'))
                        [dcmfiles, numdcm] = findFILE(dstDir, '*.img');
                    else
                        [dcmfiles, numdcm] = findFILE(dstDir, 'ref*.h5');
                    end
                    if(numdcm==0)
                        goodStatus = 0;
                    else
                        goodStatus = 1;
                    end
                    if(copy_dicom_output)
                        [dcmfiles, numdcm] = findFILE(dicomFolder, '*.dcm');
                        if(numdcm==0) 
                            goodStatus = 0;
                        end
                    end
                end
            else
                dicomFolder = fullfile(resDir, study_dates, [name '_dicom']);
                if(~isFileExist(dicomFolder))
                    goodStatus = 0;
                else
                    if(max(aif_cin_Gd_baseline_corrected(:))<1.2)
                        [dcmfiles, numdcm] = findFILE(dicomFolder, 'Image*.dcm');
                    else
                        [dcmfiles, numdcm] = findFILE(dicomFolder, 'Image_MOCO_Flow_Map_SLC0*.dcm');
                    end
                    
                    if(numdcm==0)
                        goodStatus = 0;
                    end
                end

                if goodStatus == 0
                    try
                        dstDir = fullfile(resDir, study_dates, name);
                        [resfiles, numh5] = findFILE(dstDir, 'ref*.h5');
                    catch
                        numh5 = 0;
                    end

                    if(numh5>0)
                        goodStatus = 1;
                    end
                end
            end
            disp(['Load dcm to check : ' num2str(toc)]);
            
            if(goodStatus)
                disp([num2str(n) ' out of ' num2str(num) ' - Already Processed : ' name]);
                continue;
            end
        end
        disp(['Check processed : ' num2str(toc(tstart))]);
    end
    
    try
        finfo = dir(dataName);
    catch
        disp(dataName)
        continue;
    end
    
    if(isempty(strfind(name, 'Perfusion'))~=1)
        if(finfo.bytes<100*1024*1024)
            disp(['File size too small - ' num2str(n) ' - ' name]);
            continue;
        end
    else
        try
            if(finfo.bytes<10*1024*1024)
                disp(['File size too small - ' num2str(n) ' - ' name]);
                continue;
            end
        catch
            continue;
        end
    end
    
    % start gadgetron        
    if ~exist(dstDir)
        mkdir(dstDir);
    end
    
    ts = tic;
    noise_mear_id = findNoiseDependencyMeasurementID_SavedIsmrmrd(dataName);
    disp(['find noise dependency id : ' num2str(toc(ts))]);
    
    if(~isempty(noise_mear_id))    
        noise_processed = 0;
        for kk=1:numel(noise_id_processed)
            if(strcmp(noise_id_processed{kk}, noise_mear_id)==1)
                noise_processed = 1;
                break;
            end
        end
        
        if(~noise_processed)
            ts = tic;
            [names_noise, numNoise] = findFILE(fullfile(dataDir, study_dates), ['*' noise_mear_id '*.h5']);
            disp(['find noisgit ste dependency data : ' num2str(toc(ts))]);

            if(numNoise>0)
                ts = tic;
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

                    command = ['gtdependencyquery -h ' gt_host ' -p ' GT_PORT ' -o ' fullfile(dstDir, 'noise.txt')]
                    tic; dos(command, '-echo'); disp(['nosie query : ' num2str(toc)]);
                    
                    fid = fopen(fullfile(dstDir, 'noise.txt'), 'r');
                    if(fid~=-1)
                        noise_ids = fread(fid);
                        noise_ids = char(noise_ids');
                        fclose(fid);
                        
                        if(~isempty(strfind(noise_ids, noise_mear_id)))
                            noise_processed = [noise_processed; {h5Name}];
                            noise_id_processed = [noise_id_processed; {noise_mear_id}];
                            break;
                        end
                    end
                    
                    command = ['gadgetron_ismrmrd_client -f ' h5Name ' -c default_measurement_dependencies.xml -a ' gt_host ' -p ' GT_PORT]
                    dos(command, '-echo');

                    noise_processed = [noise_processed; {h5Name}];
                    noise_id_processed = [noise_id_processed; {noise_mear_id}];
                end
                disp(['run noise dependency data : ' num2str(toc(ts))]);
            end
        else
            disp(['Noise already processed - ' noise_mear_id]);
        end
    end
            
    cd(dstDir)

    if(delete_old_res)
        tic
        delete(fullfile(dstDir, '*.h5'));
        delete(fullfile(dstDir, '*.hdr'));
        delete(fullfile(dstDir, '*.img'));
        timeUsed = toc;
        disp(['////////////////////////////////////////////////////////////////'])
        disp(['  Time to delete debugFolder on gadgetron computer (before start): ', num2str(timeUsed)])
        disp(['////////////////////////////////////////////////////////////////'])  

        if exist([dstDir,'_dicom']);
            disp('delete dicoms before processing:')
            tic; rmdir([dstDir,'_dicom'], 's'); timeUsed = toc;
            disp(['////////////////////////////////////////////////////////////////'])
            disp(['  Time to delete dicomFolder on gadgetron computer (before start): ', num2str(timeUsed)])
            disp(['////////////////////////////////////////////////////////////////'])                     
        end
    end
    
    %% run the scan
    
    files_processed = [files_processed; {name}];
    
    configNameUsed = fullfile(GTConfigFolder, configName);
   
    lenC = length(configName);
    currPath = pwd;

    lenUsed = 200 - (1+length(currPath));   
    if(lenUsed>lenC)
        lenUsed = lenC;
    end

    configNameShortened = configName(1:lenUsed);

    if( ~is_remote_computer && isunix()==0)
        debugFolder = 'D:\gtuser\mrprogs\install\DebugOutput';
        try
            mkdir(debugFolder);
        catch
        end
    else
        if(is_remote_computer)
            debugFolder = '~/Debug/DebugOutput';
            try
                mkdir(debugFolder);
            catch
            end
        else
            if(isunix())
                debugFolder = '~/Debug/DebugOutput';
            end
        end
    end
    
    if(~isempty(pre_set_debug_folder))
        debugFolder = pre_set_debug_folder;
    end
    
    if(copy_debug_output | isPerf)
        if(~is_remote_computer)
            if(isunix())
                command = ['rm -rf ' debugFolder '/*'];
                dos(command, '-echo');
            end
        else
            command = make_unix_command(user, gt_host, ['rm -rf ' debugFolder '/*.*']);
            dos(command, '-echo');
        end
    end
    
    % run the data 
    command = ['gadgetron_ismrmrd_client -f ' dataName ' -C ' configNameUsed ' -a ' gt_host ' -p ' GT_PORT ' -F ' output_format ' -G ' 'results' ' -o ref_' date_suffix '.h5']
    tic; dos(command); timeUsed = toc;
           
    if(copy_debug_output | isPerf)
        if(~is_remote_computer)
            ts = tic;
            if ~exist(dstDir); mkdir(dstDir); end
            disp(dstDir);
            if(isunix())
                ind = find(resDir=='/');
                if(isempty(ind))
                    ind = find(resDir=='/');
                end
                if(isempty(resDir_local))
                    mkdir(fullfile(dstDir, 'DebugOutput'));
                    gt_command = ['rsync -a ' debugFolder '/*.*  ' dstDir '/DebugOutput/'];
                    tic; dos(gt_command, '-echo'); timeUsed = toc;
                    disp(['////////////////////////////////////////////////////////////////'])
                    disp(['  Time to copying debugFolder to Results folder: ', num2str(timeUsed)])
                    disp(['////////////////////////////////////////////////////////////////'])

                    % fast method for deleting directory with large number of files
                    if ~exist('~/Debug/emptydir'); mkdir('~/Debug/emptydir'); end
                    command = ['rsync -a --delete ~/Debug/emptydir/ ' debugFolder,'/'];     
                    tic; dos(command, '-echo'); timeUsed = toc;
                    disp(['////////////////////////////////////////////////////////////////'])
                    disp(['  Time to delete debugFolder on gadgetron computer: ', num2str(timeUsed)])
                    disp(['////////////////////////////////////////////////////////////////'])
                else
                    %dst_dir = [resDir_local '/' resDir(ind(end)+1:end) '/' study_dates '/' name '/DebugOutput'];
                    dst_dir = [resDir_local '/' study_dates '/' name '/DebugOutput'];
                    command = ['mkdir -p ' dst_dir]
                    dos(command, '-echo');
                    command = ['mv -f ' debugFolder '/*.hdr ' dst_dir]
                    tic; dos(command, '-echo'); timeUsed = toc;
                    command = ['mv -f ' debugFolder '/*.img ' dst_dir]
                    tic; dos(command, '-echo'); timeUsed = toc;
                    disp(['////////////////////////////////////////////////////////////////'])
                    disp(['  Time to move debugFolder to Results folder: ', num2str(timeUsed)])
                    disp(['////////////////////////////////////////////////////////////////'])
                end
                
            else
                try
                    copyfile(fullfile(debugFolder, '*.*'), fullfile(dstDir, 'DebugOutput')); 
                catch
                end
                command = ['rd /s /q ' debugFolder];
                dos(command, '-echo');
            end
                        
            disp(['copy debug output : ' num2str(toc(ts))]);
        else
            ts = tic;
            ind = find(resDir=='\');
            if(isempty(ind))
                ind = find(resDir=='\');
            end
            if(isempty(resDir_local))
                dst_dir = ['/mnt/Lab-Kellman/ReconResults/' resDir(ind(end)+1:end) '/' study_dates '/' name '/DebugOutput'];
            else
                dst_dir = [resDir_local '/' resDir(ind(end)+1:end) '/' study_dates '/' name '/DebugOutput'];
            end
            
            CopyGadgetronDebugOutputOnRemote(gt_host, debugFolder, dst_dir, 0)            
            disp(['copy debug output : ' num2str(toc(ts))]);
        end
    end
    
    if(copy_dicom_output)    
        ts = tic;
        try
            [tDicom, remoteFolder] = PerformGadgetronRecon_SavedIsmrmrd_CopyDicom(resDir, name, gt_host);

            if(isempty(strfind(name, 'Perfusion'))~=1)
                missing_cases = PerformGadgetronRecon_CopyMapDicom_PerfusionCase([], {name}, resDir, []);
            end
        catch
        end
        disp(['copy dicom output : ' num2str(toc(ts))]);
    end    
    
    disp(['result folder is ' dstDir]);
    
    tUsed = [tUsed; {name, timeUsed, configName}];
end
