
function [tUsed, ignored, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, filename, gt_host, resDir, ... 
    checkProcessed, delete_old_res, startRemoteGT, configName_preset, noise_dat_processed, gt_port, copy_debug_output, copy_dicom_output, pre_set_debug_folder)
% [tUsed, ignored, noise_dat_processed] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, filename, gt_host, resDir, checkProcessed, delete_old_res, startRemoteGT, configName_preset, noise_dat_processed, gt_port, copy_debug_output, copy_dicom_output, pre_set_debug_folder)
% [tUsed, ignored] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData('I:\KAROLINSKA', 'xxxx', 'localhost', 'I:\ReconResults\KAROLINSKA')
% setenv('OutputFormat', 'h5')

if(strcmp(gt_host, 'gt1'))
    gt_host = '137.187.135.169';
end
if(strcmp(gt_host, 'gt2'))
    gt_host = '137.187.135.238';
end
if(strcmp(gt_host, 'gt3'))
    gt_host = '137.187.135.155';
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
    
%     if(copy_debug_output)
%         isPerf = 1;
%     end
       
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
    
%     if(str2num(measurementID)>300000)
%         continue;
%     end
    
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
        if(finfo.bytes<120*1024*1024)
            disp(['File size too small - ' num2str(n) ' - ' name]);
%             cd(dataDir)
            ignored = [ignored; {n, name, finfo.bytes/1024}];
            continue;
        end
    else
        try
            if(finfo.bytes<5*1024*1024)
                disp(['File size too small - ' num2str(n) ' - ' name]);
    %             cd(dataDir)
                ignored = [ignored; {n, name, finfo.bytes/1024}];
                continue;
            end
        catch
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
                %% 
            end
        end
        disp(['Start remote gadgetron : ' num2str(toc(tstart))]);
    end
    
    
    if ~exist(dstDir); mkdir(dstDir); end    
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
        if ~exist('~/Debug/emptydir'); mkdir('~/Debug/emptydir'); end
        disp('delete DebugOutput folder before processing:')
        command = ['sudo rsync -a --delete ~/Debug/emptydir/ ' dstDir,'/DebugOutput/'];     
        dos(command, '-echo');   
        command = ['sudo rsync -a --delete ~/Debug/emptydir/ ' dstDir,'/'];     
        %dos(command, '-echo');   
        tic; dos(command, '-echo'); timeUsed = toc;
        disp(['////////////////////////////////////////////////////////////////'])
        disp(['  Time to delete debugFolder on gadgetron computer (before start): ', num2str(timeUsed)])
        disp(['////////////////////////////////////////////////////////////////'])  
        % delete dicoms
        if exist([dstDir,'_dicom']);
            disp('delete dicoms before processing:')
            command = ['sudo rsync -a --delete ~/Debug/emptydir/ ' dstDir,'_dicom/'];     
            %dos(command, '-echo'); 
            tic; dos(command, '-echo'); timeUsed = toc;
            disp(['////////////////////////////////////////////////////////////////'])
            disp(['  Time to delete dicomFolder on gadgetron computer (before start): ', num2str(timeUsed)])
            disp(['////////////////////////////////////////////////////////////////'])         
            
        end
        
        if 0 % remove old code
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
        %delete(fullfile(dstDir,filesep,'DebugOutput',filesep, '*'));  
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
        try
            if ~exist(debugFolder); mkdir(debugFolder); end
        catch
        end
    end
    
    if(copy_debug_output | isPerf)
        if(~is_remote_computer)
            if(isunix())
                command = ['sudo rm -rf ' debugFolder '/*.*'];
                dos(command, '-echo');
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
           
    if(copy_debug_output | isPerf)
        if(~is_remote_computer)
            ts = tic;
            if ~exist(dstDir); mkdir(dstDir); end
            disp(dstDir);
            if(isunix())
                mkdir(fullfile(dstDir, 'DebugOutput'));
                %gt_command = ['cp -r ' fullfile(debugFolder, '*.*') ' ' dstDir '/DebugOutput'];
                gt_command = ['rsync -a ' debugFolder '/*.*  ' dstDir '/DebugOutput/'];
                %dos(gt_command, '-echo');
                tic; dos(gt_command, '-echo'); timeUsed = toc;
                disp(['////////////////////////////////////////////////////////////////'])
                disp(['  Time to copying debugFolder to Results folder: ', num2str(timeUsed)])
                disp(['////////////////////////////////////////////////////////////////'])
                
                %command = ['sudo rm -rf ' debugFolder '/*.*'];
                %dos(command, '-echo');
                
                % fast method for deleting directory with large number of files
                if ~exist('~/Debug/emptydir'); mkdir('~/Debug/emptydir'); end
                command = ['sudo rsync -a --delete ~/Debug/emptydir/ ' debugFolder,'/'];     
                tic; dos(command, '-echo'); timeUsed = toc;
                disp(['////////////////////////////////////////////////////////////////'])
                disp(['  Time to delete debugFolder on gadgetron computer: ', num2str(timeUsed)])
                disp(['////////////////////////////////////////////////////////////////'])
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
            %debugFolder = ['/home/' user '/Debug/DebugOutput']
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
            CopyGadgetronDebugOutputOnRemote(gt_host, debugFolder, dst_dir, 1)
            
            disp(['copy debug output : ' num2str(toc(ts))]);
        end
    end
    
    if(copy_dicom_output)    
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

            if(startRemoteGT)
                copyfile(['D:\Temp\record_' GT_PORT '.txt'], [dstDir '\record_' gt_host '_' GT_PORT '.txt']);
            end
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
    
    disp(['result folder is ' dstDir]);
    
    tUsed = [tUsed; {name, timeUsed, configName}];
end
