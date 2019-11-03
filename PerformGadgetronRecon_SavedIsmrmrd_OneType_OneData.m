
function [tUsed, ignored, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, filename, gt_host, resDir, ... 
    checkProcessed, delete_old_res, startRemoteGT, configName_preset, noise_dat_processed, gt_port, copy_debug_output)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, filename, gt_host, resDir, checkProcessed, delete_old_res, startRemoteGT, configName_preset, noise_dat_processed)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData('I:\KAROLINSKA', 'xxxx', 'localhost', 'I:\ReconResults\KAROLINSKA')
% setenv('OutputFormat', 'h5')

if(strcmp(gt_host, 'gt1'))
    gt_host = '137.187.135.97';
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
    checkProcessed = 1;
end

if(nargin<6)
    delete_old_res = 1;
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
    gt_port = GT_PORT;
end

if(nargin<11)
    copy_debug_output = 0;
end

setenv('GT_PORT', gt_port);
GT_PORT = gt_port;

GTHome = getenv('GADGETRON_HOME');
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

xmlUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xml';

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
    
    if(copy_debug_output)
        isPerf = 1;
    end
       
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
    
    dataName = fullfile(dataDir, study_dates, [name '.h5']);

%     if(strcmp(gt_host, 'localhost')==1 && isunix()==0)
%         [pathstr, configName, ext] = fileparts(configName);
%         configName = [configName '_localhost' ext];
%     end
    
    if(~isempty(configName_preset{1}))
        if(n>numel(configName_preset))
            configName = configName_preset{end};
        else
            configName = configName_preset{n};
        end
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
                    % a = readGTPlusExportImageSeries_Squeeze(dstDir, 120);
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
%                     
%                     if(max(aif_cin_Gd_baseline_corrected)>1.5)
%                         dstAcqTimes = analyze75read(fullfile(dstDir, 'DebugOutput', 'dstAcqTimes_0.hdr'));
%                         AIF_AcqTimes = analyze75read(fullfile(dstDir, 'DebugOutput', 'AIF_AcqTimes_0.hdr'));
%                         perf_LUT = analyze75read(fullfile(dstDir, 'DebugOutput', 'Perf_T1_Correction_LUT.hdr'));
% 
%                         sr_0 = analyze75read(fullfile(dstDir, 'DebugOutput', 'input_for_SRNorm_0.hdr'));
%                         pd_0 = analyze75read(fullfile(dstDir, 'DebugOutput', 'inputPD_for_SRNorm_0.hdr'));
%                         sr_norm_0 = analyze75read(fullfile(dstDir, 'DebugOutput', 'SRNorm_0.hdr'));
%                         pd0 = analyze75read(fullfile(dstDir, 'DebugOutput', 'PD_0.hdr'));
%                         perf_mask_0 = analyze75read(fullfile(dstDir, 'DebugOutput', 'perf_mask_0.hdr'));
%                         gd0 = analyze75read(fullfile(dstDir, 'DebugOutput', 'CASignal_Perf_0.hdr'));
%                         Perf_AcqTimes_0 = analyze75read(fullfile(dstDir, 'DebugOutput', 'Perf_AcqTimes_0.hdr'));
%                         gd0_resampled = analyze75read(fullfile(dstDir, 'DebugOutput', 'Input_perf_computeFlowMap_0.hdr'));
%                     end

%                     sr_1 = analyze75read(fullfile(dstDir, 'DebugOutput', 'input_for_SRNorm_1.hdr'));
%                     pd_1 = analyze75read(fullfile(dstDir, 'DebugOutput', 'inputPD_for_SRNorm_1.hdr'));
%                     sr_norm_1 = analyze75read(fullfile(dstDir, 'DebugOutput', 'SRNorm_1.hdr'));
%                     pd1 = analyze75read(fullfile(dstDir, 'DebugOutput', 'PD_1.hdr'));
%                     perf_mask_1 = analyze75read(fullfile(dstDir, 'DebugOutput', 'perf_mask_1.hdr'));
%                     gd1 = analyze75read(fullfile(dstDir, 'DebugOutput', 'CASignal_Perf_1.hdr'));
%                     Perf_AcqTimes_1 = analyze75read(fullfile(dstDir, 'DebugOutput', 'Perf_AcqTimes_1.hdr'));
%                     gd1_resampled = analyze75read(fullfile(dstDir, 'DebugOutput', 'Input_perf_computeFlowMap_1.hdr'));
% 
%                     sr_2 = analyze75read(fullfile(dstDir, 'DebugOutput', 'input_for_SRNorm_2.hdr'));
%                     pd_2 = analyze75read(fullfile(dstDir, 'DebugOutput', 'inputPD_for_SRNorm_2.hdr'));
%                     sr_norm_2 = analyze75read(fullfile(dstDir, 'DebugOutput', 'SRNorm_2.hdr'));
%                     pd2 = analyze75read(fullfile(dstDir, 'DebugOutput', 'PD_2.hdr'));
%                     perf_mask_2 = analyze75read(fullfile(dstDir, 'DebugOutput', 'perf_mask_2.hdr'));
%                     gd2 = analyze75read(fullfile(dstDir, 'DebugOutput', 'CASignal_Perf_2.hdr'));
%                     Perf_AcqTimes_2 = analyze75read(fullfile(dstDir, 'DebugOutput', 'Perf_AcqTimes_2.hdr'));
%                     gd2_resampled = analyze75read(fullfile(dstDir, 'DebugOutput', 'Input_perf_computeFlowMap_2.hdr'));
                    disp(['Load hdr to check : ' num2str(toc)]);
                catch
                    goodStatus = 0;
                end
            end

%             [hdrfiles, numhdr] = findFILE(dstDir, '*.img');
                          
            tic
            if(~isPerf)
                hdrFolder = fullfile(resDir, study_dates, [name]);
                dicomFolder = fullfile(resDir, study_dates, [name '_dicom']);
                if(~isFileExist(dicomFolder))
                    goodStatus = 0;
                    if(isFileExist(hdrFolder))
                        try
                        [dcmfiles, numdcm] = findFILE(dstDir, '*.img');
                        catch
                            numdcm = 0;
                        end
                        if(numdcm==0)
                            goodStatus = 0;
                        else
                            goodStatus = 1;
                        end
                    end
                else
                    try
                        [dcmfiles, numdcm] = findFILE(dstDir, '*.img');
                    catch
                        numdcm = 0;
                    end
                    if(numdcm==0)
                        goodStatus = 0;
                    else
                        goodStatus = 1;
                    end
                end
            else
%                 dicomFolder = fullfile(resDir, study_dates, [name '_dicom']);
%                 if(~isFileExist(dicomFolder))
%                     goodStatus = 0;
%                 else
%                     if(max(aif_cin_Gd_baseline_corrected(:))<1.2)
%                         [dcmfiles, numdcm] = findFILE(dicomFolder, 'Image*.dcm');
%                     else
%                         [dcmfiles, numdcm] = findFILE(dicomFolder, 'Image_MOCO_Flow_Map_SLC0*.dcm');
%                     end
%                     
%                     if(numdcm==0)
%                         goodStatus = 0;
%                     end
%                 end
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
    
%     isPerf = 0;
    if(isempty(strfind(name, 'Perfusion'))~=1)
%         isPerf = 1;
        if(finfo.bytes<100*1024*1024)
            disp(['File size too small - ' num2str(n) ' - ' name]);
%             cd(dataDir)
            ignored = [ignored; {n, name, finfo.bytes/1024}];
            continue;
        end
    else
        if(finfo.bytes<10*1024*1024)
            disp(['File size too small - ' num2str(n) ' - ' name]);
%             cd(dataDir)
            ignored = [ignored; {n, name, finfo.bytes/1024}];
            continue;
        end
    end
    
    % start gadgetron
    if(startRemoteGT)
        tstart = tic;
        if(~is_remote_computer && isunix()==0)
            cd('D:\gtuser\gt_scanner_setup_scripts')
            command = ['gadgetron -p %GT_PORT% > D:\Temp\record_' GT_PORT '.txt']
            dos([command ' &'])
        else
            if (~isempty(user) & startRemoteGT)
                StopGadgetronOnRemote(gt_host, GT_PORT);                
                StartGadgetronOnRemote(gt_host, GT_PORT);
            end
        end
        disp(['Start remote gadgetron : ' num2str(toc(tstart))]);
    end
    
    mkdir(dstDir);
    
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
        disp(['delete dstDir : starts ... ']);
        ts = tic;
        delete(fullfile(dstDir, 'res*.h5'));
        delete(fullfile(dstDir, 'out*.h5'));
        delete(fullfile(dstDir, 'ref*.h5'));
        delete(fullfile(dstDir, '*.xml'));

        delete(fullfile(dstDir, '*.nii'));
%         delete(fullfile(dstDir, 'gadgetron_*.hdr'));
%         delete(fullfile(dstDir, 'gadgetron_*.img'));
%         delete(fullfile(dstDir, 'Generic*.hdr'));
%         delete(fullfile(dstDir, 'Generic*.img'));
%         delete(fullfile(dstDir, 'GTPrep*.hdr'));
%         delete(fullfile(dstDir, 'GTPrep*.img'));
%         delete(fullfile(dstDir, 'GT*.hdr'));
%         delete(fullfile(dstDir, 'GT*.img'));
        delete(fullfile(dstDir, '*.attrib'));
        delete(fullfile(dstDir, '*.hdr'));
        delete(fullfile(dstDir, '*.img'));
        
        delete(fullfile(dstDir, '*.xml'));                 
        
        dicomDir = fullfile(resDir, study_dates, [name '_dicom']);
        delete(fullfile(dicomDir, '*.*'));   
        disp(['delete dstDir : ' num2str(toc(ts))]);
    end
    
    %% run the scan
    
%     disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' name]);
    % pause(0.5)
    
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
        if( ~is_remote_computer && isunix()==0)
            debugFolder = 'D:\gtuser\mrprogs\install\DebugOutput';
%             try
%                 rmdir(debugFolder, 's');
%             catch
%             end

            try
                mkdir(debugFolder);
            catch
            end
        end
        
        if(~is_remote_computer && isunix()==1)
            debugFolder = '~/Debug/DebugOutput';
%             try
%                 rmdir(debugFolder, 's');
%             catch
%             end

            try
                mkdir(debugFolder);
            catch
            end
        end
    end
    
    % run the data 
    if(isPerf)
%         command = ['gadgetron_ismrmrd_client -f ' dataName ' -C ' configNameUsed ' -a ' gt_host ' -p ' GT_PORT ' -F ' output_format ' -G ' configNameShortened ' -o ref_' date_suffix '.h5']
        command = ['gadgetron_ismrmrd_client -f ' dataName ' -C ' configNameUsed ' -a ' gt_host ' -p ' GT_PORT ' -F ' output_format ' -G ' 'results' ' -o ref_' date_suffix '.h5']
    else
%         command = ['gadgetron_ismrmrd_client -f ' dataName ' -C ' configNameUsed ' -a ' gt_host ' -p ' GT_PORT ' -F ' output_format ' -G ' configNameShortened ' -o ref_' date_suffix '.h5']
        command = ['gadgetron_ismrmrd_client -f ' dataName ' -C ' configNameUsed ' -a ' gt_host ' -p ' GT_PORT ' -F ' output_format ' -G ' 'results' ' -o ref_' date_suffix '.h5']
    end
    tic; dos(command); timeUsed = toc;
           
    if(isPerf)
        if(~is_remote_computer)
            ts = tic;
            mkdir(dstDir);
            disp(dstDir);
            if(isunix())
                mkdir(fullfile(dstDir, 'DebugOutput'));
                gt_command = ['cp -r ' fullfile(debugFolder, '*.*') ' ' dstDir '/DebugOutput'];
                gt_command
                dos(gt_command, '-echo');
            else
                copyfile(fullfile(debugFolder, '*.*'), fullfile(dstDir, 'DebugOutput')); 
            end
            
%             command = ['move /Y ' debugFolder ' ' dstDir];
%             dos(command, '-echo');
            disp(['copy debug output : ' num2str(toc(ts))]);
        else
            debug_folder = ['/home/' user '/Debug/DebugOutput']
            ts = tic;

            try
                ind = find(resDir=='/');
                if(isempty(ind))
                    ind = find(resDir=='/');
                end
                dst_dir = ['/mnt/Lab-Kellman/ReconResults/' resDir(ind(end)+1:end) '/' study_dates '/' name '/DebugOutput'];
            catch
                ind = find(resDir=='\');
                if(isempty(ind))
                    ind = find(resDir=='\');
                end
                dst_dir = ['/mnt/Lab-Kellman/ReconResults/' resDir(ind(end)+1:end) '/' study_dates '/' name '/DebugOutput'];
            end
            try
                mkdir(dst_dir);
            catch
            end
            CopyGadgetronDebugOutputOnRemote(gt_host, debug_folder, dst_dir, 1)
            
            disp(['copy debug output : ' num2str(toc(ts))]);
        end
    end
    
    if(strcmp(gt_host, '137.187.134.169')==0)    
        ts = tic;
        try
            [tDicom, remoteFolder] = PerformGadgetronRecon_SavedIsmrmrd_CopyDicom(resDir, name, gt_host);

    %         [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
    %         dicomDir = fullfile(resDir, study_dates, [name '_dicom']);
    %         rmdir(dicomDir)

            if(isempty(strfind(name, 'Perfusion'))~=1)
                missing_cases = PerformGadgetronRecon_CopyMapDicom_PerfusionCase([], {name}, resDir, []);
            end
        catch
        end
        disp(['copy dicom output : ' num2str(toc(ts))]);

        ts = tic;
        dstDir = fullfile(resDir, study_dates, name);   
        if(~is_remote_computer && isunix()==0)
            if(startRemoteGT)
                command = ['taskkill /F /FI "IMAGENAME eq gadgetron.*"'];
                dos(command)

                command = ['taskkill /F /IM cmd.exe'];
                dos(command)
            end

            copyfile(['D:\Temp\record_' GT_PORT '.txt'], [dstDir '\record_' gt_host '_' GT_PORT '.txt']);
        else
            if (~isempty(user) & startRemoteGT)
                StopGadgetronOnRemote(gt_host, GT_PORT);
                CopyGadgetronRecordOnRemote(gt_host, GT_PORT, [dstDir '/record_' gt_host '_' GT_PORT '.txt']);
            end
        end
        disp(['shutdown gadgetron : ' num2str(toc(ts))]);

    %     if(sendDicom)
    %         dicomServer = 'barbados';
    %         dicomPort = 11112;
    %         command = ['D:\gtuser\gt_windows_setup\dcmtk-3.6.0\install_vc14\bin\storescu ' dicomServer '.nhlbi.nih.gov ' num2str(dicomPort) ' ' remoteFolder ' --scan-directories -aec DCM4CHEE --user admin --password admin'];
    %         tic; dos(command, '-echo'); toc
    %     end
    end    
    
    tUsed = [tUsed; {name, timeUsed, configName}];
end
