function timeUsed = run_gt_recon(folderDir, dataName, h5Name, deleteh5, isVD, isVD11, isNX, isNX20, isAdjScan, configName, resDir, styleSheet, startRemoteGT, h5Only, remoteXml, compressionBit, paraXml, debug_folder)
% timeUsed = run_gt_recon(folderDir, dataName, h5Name, deleteh5, isVD, isVD11, isNX, isNX20, isAdjScan, configName, resDir, styleSheet, startRemoteGT, h5Only, remoteXml, compressionBit, paraXml, debug_folder)

% siemens to hdf5

currDir = pwd;

GTHome = getenv('GADGETRON_HOME')
%GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
if(~isunix())
    GTConfigFolder = fullfile([GTHome '_master_debug'], 'share/gadgetron/config');
else
    GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
end
disp(['GTConfigFolder is ' GTConfigFolder]);
if(remoteXml)
    configNameUsed = configName;
else
    configNameUsed = fullfile(GTConfigFolder, configName);
end
disp(['----> configNameUsed is ' configNameUsed]);

GT_HOST = getenv('GT_HOST')
GT_PORT = getenv('GT_PORT')
OutputFormat = getenv('OutputFormat')

if(startRemoteGT)
    if((strcmp(GT_HOST, 'localhost')==1))
        if(isunix())
            command = ['gadgetron -p ' GT_PORT ' > ~/record_' GT_PORT '.txt']
            dos([command ' &'])
        else            
            cd('D:/gtuser/gt_scanner_setup_scripts')
            command = ['gadgetron -p ' GT_PORT ' > D:/Temp/record_' GT_PORT '.txt']
            dos([command ' &'])
        end
    else
        [key, user] = sshKeyLookup(GT_HOST);
        if (~isempty(user) & startRemoteGT)
            StartGadgetronOnRemote(GT_HOST, GT_PORT);
        end
    end
end

if ( ~exist(resDir) )
    mkdir(resDir);
end
cd(resDir)

if(~h5Only)
    if ( ~isFileExist(h5Name) || deleteh5 )
        delete(h5Name);
    %     command = ['siemens_to_HDF5 ' dataName ' ' h5Name]
    %     dos(command);
    end
end

if(isunix())
    schema_dir = [GTHome '/schema/'];
else
    schema_dir = '%GADGETRON_DIR%/install/schema';
end

if(isNX20)
    styleSheetUsed = [schema_dir '/wip_071_qPerf_IsmrmrdParameterMap_Siemens_Perfusion_NX20.xsl'];
end

if(isNX)
    %styleSheetUsed = [schema_dir '/wip_070_fire_IsmrmrdParameterMap_Siemens.xsl'];
    styleSheetUsed = [schema_dir '/IsmrmrdParameterMap_Siemens.xsl'];
else
    styleSheetUsed = [schema_dir '/IsmrmrdParameterMap_Siemens.xsl'];
end

if ( nargin >= 8 )
    styleSheetUsed = [schema_dir '/' styleSheet];
end

flashRef = 0;
if( ~isempty(strfind(styleSheet, 'IsmrmrdParameterMap_Siemens_EPI_FLASHREF')) )
    flashRef = 1;
end

if(isNX20)
    xmlUsed = [schema_dir '/wip_071_qPerf_IsmrmrdParameterMap_Siemens_Perfusion_NX20.xml'];
elseif(isNX)
    %xmlUsed = [schema_dir '/wip_070_fire_IsmrmrdParameterMap_Siemens.xml'];
    xmlUsed = [schema_dir '/IsmrmrdParameterMap_Siemens.xml'];
else
    xmlUsed = [schema_dir '/IsmrmrdParameterMap_Siemens_Perfusion.xml'];
end

if(~isempty(paraXml))
    if(nargin>=16)
        xmlUsed = [schema_dir '/' paraXml];
    end
end

% dependency_xml = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml';
dependency_xml = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml';

% if VD, run the dependent scan

lenC = length(configName);
currPath = pwd;

lenUsed = 200 - (1+length(currPath));   
if(lenUsed>lenC)
    lenUsed = lenC;
end

configNameShortened = configName(1:lenUsed);

date_suffix = datestr(date, 'yyyymmdd');
time_now = datestr(now, 'HHMMSS');
date_suffix = [date_suffix '_' time_now];

xml_opt = ' -C ';
if(remoteXml)
    xml_opt = ' -c ';
end

compression_opt = [];
if (compressionBit>0)
    if(compressionBit>1)
        compression_opt = [' -P ' num2str(compressionBit) ' '];
    else
        compression_opt = [' -T ' num2str(compressionBit) ' '];
    end
end

studyDate = datestr(date, 'yyyy-mm-dd');

if ( isVD | isNX | isNX20 )
    if ( isVD11 )
        if ( ~isFileExist(h5Name) || deleteh5 )
            delete(h5Name);
            command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 -X -H --studyDate ' studyDate]
            dos(command);
            
            command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 --studyDate ' studyDate]
            dos(command);
        end
        
        command = ['gadgetron_ismrmrd_client -f ' h5Name xml_opt configNameUsed ' -a ' GT_HOST ' -p ' GT_PORT ' -F ' OutputFormat ' -G ' configNameShortened ' -o ref_' date_suffix '.h5' compression_opt]
        command
        tic; dos(command); timeUsed = toc;
        
        if(~isempty(debug_folder))
            mkdir(fullfile(resDir, 'DebugOutput'));
            movefile(fullfile(debug_folder, '*.*'), fullfile(resDir, 'DebugOutput'));
        end
    else

        if(h5Only)
            
            [names, num] = findFILE(folderDir, 'ISMRMRD_Noise_dependency_*.h5');   
            if(num==0)
                [names, num] = findFILE(folderDir, '*noise*.h5');   
            end
            for n=1:num
                command = ['gadgetron_ismrmrd_client -f ' names{n} ' -c ' dependency_xml ' -a ' GT_HOST ' -p ' GT_PORT]
                dos(command, '-echo');
            end            
            
            command = ['gadgetron_ismrmrd_client -f ' h5Name xml_opt configNameUsed ' -a ' GT_HOST ' -p ' GT_PORT ' -F ' OutputFormat ' -G ' configNameShortened ' -o ref_' date_suffix '.h5' compression_opt]                           
            command
            tic; dos(command); timeUsed = toc;
            
            if(~isempty(debug_folder))
                try
                    mkdir(fullfile(resDir, 'DebugOutput'));
                    movefile(fullfile(debug_folder, '*.*'), fullfile(resDir, 'DebugOutput'));
                catch
                end
            end
        else
            hasAdj = 1;
            if(deleteh5)
                % check whether two measurements are in the dat file
                command = ['siemens_to_ismrmrd -f  ' dataName ' -o test_meas.h5 --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 --studyDate ' studyDate]
                dos(command, '-echo');

                if(~isFileExist('test_meas.h5'))
                    hasAdj = 0;
                else
                    delete('test_meas.h5');
                end
            end

            if(isNX || isNX20)
                noise_xsl = 'IsmrmrdParameterMap_Siemens_NX.xsl';
            else
                noise_xsl = 'IsmrmrdParameterMap_Siemens.xsl';
            end
            
            [names, num] = findFILE(folderDir, '*AdjCoilSens*.dat');          
            for n=1:num
                command = ['siemens_to_ismrmrd -f ' names{n} ' -o ' names{n} '_AdjCoilSens.h5 --user-map ' schema_dir '/IsmrmrdParameterMap_Siemens.xml --user-stylesheet ' schema_dir '/' noise_xsl ' -z 1 --studyDate ' studyDate]
                dos(command, '-echo');
                
                command = ['gadgetron_ismrmrd_client -f ' names{n} '_AdjCoilSens.h5 -c ' dependency_xml ' -a ' GT_HOST ' -p ' GT_PORT]
                dos(command, '-echo');
            end            
            
            has_only_one_scan = 0;
            delete('config_buffer*.*');
            delete('*.xml');
            command = ['siemens_to_ismrmrd -f  ' dataName ' -o test.h5 --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 -X -H --studyDate ' studyDate];
            dos(command, '-echo');
            delete('test.h5');
            if(~exist('config_buffer.xprot'))
                has_only_one_scan = 1;
                disp(['This data has only one scan !!! ']);
                hasAdj = 0;
            end
            
            if(hasAdj)
                if(~has_only_one_scan)
                    if ( ~isFileExist('AdjCoilSens.h5') || deleteh5 )
                        delete('AdjCoilSens.h5');
                        delete('*.*');
                        command = ['siemens_to_ismrmrd -f ' dataName ' -o AdjCoilSens.h5 --skipSyncData --user-map ' schema_dir '/IsmrmrdParameterMap_Siemens_NX.xml --user-stylesheet ' schema_dir '/' noise_xsl ' -z 1 --studyDate ' studyDate]
                        dos(command, '-echo');
                    end

                    command = ['gadgetron_ismrmrd_client -f AdjCoilSens.h5  -c ' dependency_xml ' -a ' GT_HOST ' -p ' GT_PORT]
                    dos(command, '-echo');
                end
            end

            if ( ~isFileExist(h5Name) || deleteh5 )

                delete(h5Name);

                if(hasAdj)       
                    command = ['siemens_to_ismrmrd -f  ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 3 -X -H --studyDate ' studyDate]
                    if(flashRef)
                        command = [command ' -F'];
                    end
                    if(isNX)
                        command = [command ' --skipSyncData'];
                    end
                    dos(command, '-echo');

                    command = ['siemens_to_ismrmrd -f  ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 3 --studyDate ' studyDate]
                    if(flashRef)
                        command = [command ' -F'];
                    end
                    if(isNX)
                        command = [command ' --skipSyncData'];
                    end
                    dos(command, '-echo');
                    
                    if(~isFileExist(h5Name))                    
                        command = ['siemens_to_ismrmrd -f  ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 -X -H --studyDate ' studyDate]
                        if(flashRef)
                            command = [command ' -F'];
                        end
                        if(isNX)
                            command = [command ' --skipSyncData'];
                        end
                        dos(command, '-echo');

                        command = ['siemens_to_ismrmrd -f  ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 --studyDate ' studyDate]
                        if(flashRef)
                            command = [command ' -F'];
                        end
                        if(isNX)
                            command = [command ' --skipSyncData'];
                        end
                        dos(command, '-echo');
                    end                    
                else
                    command = ['siemens_to_ismrmrd -f  ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 1 -X -H --studyDate ' studyDate]
                    if(flashRef)
                        command = [command ' -F'];
                    end
                    if(isNX)
                        command = [command ' --skipSyncData'];
                    end
                    dos(command, '-echo');

                    command = ['siemens_to_ismrmrd -f  ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 1 --studyDate ' studyDate]
                    if(flashRef)
                        command = [command ' -F'];
                    end
                    if(isNX)
                        command = [command ' --skipSyncData'];
                    end
                    dos(command, '-echo');
                end
            end

            command = ['gadgetron_ismrmrd_client -f ' h5Name xml_opt configNameUsed ' -a ' GT_HOST ' -p ' GT_PORT ' -F ' OutputFormat ' -G ' configNameShortened ' -o ref_' date_suffix '.h5' compression_opt]
            command
            tic; dos(command); timeUsed = toc;
            
            if(~isempty(debug_folder))
                mkdir(fullfile(resDir, 'DebugOutput'));
                try
                    movefile(fullfile(debug_folder, '*.*'), fullfile(resDir, 'DebugOutput'));
                catch
                end
            end
        end
    end
else  
    if(isAdjScan)
        command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' schema_dir '/IsmrmrdParameterMap_Siemens.xml --user-stylesheet ' schema_dir '/IsmrmrdParameterMap_Siemens.xsl -z 1 -X -H --studyDate ' studyDate]
        dos(command);
            
        command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' schema_dir '/IsmrmrdParameterMap_Siemens.xml --user-stylesheet ' schema_dir '/IsmrmrdParameterMap_Siemens.xsl -z 1 -X --studyDate ' studyDate]
        dos(command);
            
        command = ['gadgetron_ismrmrd_client -f ' h5Name ' -c ' dependency_xml ' -a ' GT_HOST ' -p ' GT_PORT]
        tic; dos(command); timeUsed = toc;
    else
        if ( ~isFileExist(h5Name) || deleteh5 )
            delete(h5Name);
            command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' schema_dir '/IsmrmrdParameterMap_Siemens_VB17.xml --user-stylesheet ' styleSheetUsed ' -z 1 -X -H --studyDate ' studyDate]
            dos(command);

            command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' schema_dir '/IsmrmrdParameterMap_Siemens_VB17.xml --user-stylesheet ' styleSheetUsed ' -z 1 --studyDate ' studyDate]
            dos(command);
        end

        if ~exist(configNameUsed)
            xml_opt = ' -c ';
            [xml_path, xml_name, ext] = fileparts(configNameUsed);
            configNameUsed = [xml_name ext];
        end

        command = ['gadgetron_ismrmrd_client -f ' h5Name xml_opt configNameUsed ' -a ' GT_HOST ' -p ' GT_PORT ' -F ' OutputFormat ' -G ' configNameShortened ' -o ref_' date_suffix '.h5' compression_opt]
        command
        tic; dos(command); timeUsed = toc;
        
        if(~isempty(debug_folder))
            mkdir(fullfile(resDir, 'DebugOutput'));
            movefile(fullfile(debug_folder, '*.*'), fullfile(resDir, 'DebugOutput'));
        end
    end
end

if(strcmp(GT_HOST, 'localhost')==1)
    if(startRemoteGT)
        if(isunix())
            StopGadgetronOnRemote(GT_HOST, GT_PORT);
        else
            command = ['taskkill /F /FI "IMAGENAME eq gadgetron.*"'];
            dos(command)
        end
    end
    
    try
        copyfile(['D:/Temp/record_' GT_PORT '.txt'], [resDir '/record_' GT_HOST '_' GT_PORT '.txt']);
    catch
    end
else
    [key, user] = sshKeyLookup(GT_HOST);
    if (~isempty(user) & startRemoteGT)
        StopGadgetronOnRemote(GT_HOST, GT_PORT);
        CopyGadgetronRecordOnRemote(GT_HOST, GT_PORT, [resDir '/record_' GT_HOST '_' GT_PORT '.txt']);
    end
end

cd(currDir)

end

function res = PerformGadgetronRecon_Matlab_ROIValues_OneCase(s1, s2, s3, a1, a2, a3)
%  res = PerformGadgetronRecon_Matlab_ROIValues_OneCase(s1, s2, s3, a1, a1, a3)
% given the ROIs s1/s2/s3, get the flow and other values
% res has [flow, Ki, E, Visf, Vp, PS, SD, Ki_MF, Ki_Fermi, Ki_TwoCompExp, Ki_BTEX]

if(isfield(a1, 'flowmaps_grappa_PSIR'))
    %% stress
       
    [f1, f2, f3] = get_roi_values(a1.flowmaps_grappa_PSIR, a2.flowmaps_grappa_PSIR, a3.flowmaps_grappa_PSIR, s1, s2, s3);
    res.flow = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(stress.E_stress, s1, s2, s3);
    res.E = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(stress.Visf_stress, s1, s2, s3);
    res.Visf = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(stress.Vp_stress, s1, s2, s3);
    res.Vp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(stress.PS_stress, s1, s2, s3);
    res.PS = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(stress.SDMap_stress, s1, s2, s3);
    res.SD = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(stress.Ki_stress(:,:,1,:)), s1, s2, s3);
    res.Ki_MF = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(stress.Ki_stress(:,:,2,:)), s1, s2, s3);
    res.Ki_Fermi = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(stress.Ki_stress(:,:,3,:)), s1, s2, s3);
    res.Ki_TwoCompExp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(stress.Ki_stress(:,:,4,:)), s1, s2, s3);
    res.Ki_BTEX = [f1.m f2.m f3.m];

    %% if having second roi
    if(numel(s1.ROI_info_table)==2)

        [f1, f2, f3] = get_2nd_roi_values(stress.flow_stress, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.E_stress, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.Visf_stress, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.Vp_stress, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.PS_stress, s1, s2, s3);
        res.PS_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.SDMap_stress, s1, s2, s3);
        res.SD_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(stress.Ki_stress(:,:,1,:)), s1, s2, s3);
        res.Ki_MF_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(stress.Ki_stress(:,:,2,:)), s1, s2, s3);
        res.Ki_Fermi_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(stress.Ki_stress(:,:,3,:)), s1, s2, s3);
        res.Ki_TwoCompExp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(stress.Ki_stress(:,:,4,:)), s1, s2, s3);
        res.Ki_BTEX_i = [f1.m f2.m f3.m];
    end
else
    %% rest
    
    rest = a;
    
    [f1, f2, f3] = get_roi_values(rest.flow_rest, s1, s2, s3);
    res.flow = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(rest.E_rest, s1, s2, s3);
    res.E = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(rest.Visf_rest, s1, s2, s3);
    res.Visf = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(rest.Vp_rest, s1, s2, s3);
    res.Vp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(rest.PS_rest, s1, s2, s3);
    res.PS = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(rest.SDMap_rest, s1, s2, s3);
    res.SD = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(rest.Ki_rest(:,:,1,:)), s1, s2, s3);
    res.Ki_MF = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(rest.Ki_rest(:,:,2,:)), s1, s2, s3);
    res.Ki_Fermi = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(rest.Ki_rest(:,:,3,:)), s1, s2, s3);
    res.Ki_TwoCompExp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(rest.Ki_rest(:,:,4,:)), s1, s2, s3);
    res.Ki_BTEX = [f1.m f2.m f3.m];

    %% if having second roi
%     if(numel(s1.ROI_info_table)==2)

        [f1, f2, f3] = get_2nd_roi_values(rest.flow_rest, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(rest.E_rest, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(rest.Visf_rest, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(rest.Vp_rest, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(rest.PS_rest, s1, s2, s3);
        res.PS_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(rest.SDMap_rest, s1, s2, s3);
        res.SD_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(rest.Ki_rest(:,:,1,:)), s1, s2, s3);
        res.Ki_MF_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(rest.Ki_rest(:,:,2,:)), s1, s2, s3);
        res.Ki_Fermi_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(rest.Ki_rest(:,:,3,:)), s1, s2, s3);
        res.Ki_TwoCompExp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(rest.Ki_rest(:,:,4,:)), s1, s2, s3);
        res.Ki_BTEX_i = [f1.m f2.m f3.m];
%     end
end

end

function [f1, f2, f3] = get_roi_values(a1, a2, a3, s1, s2, s3)
    f1 = roi_statistics(a1(:,:,end), s1.ROI_info_table(1,1));
    f2 = roi_statistics(a2(:,:,end), s2.ROI_info_table(1,1));
    f3 = roi_statistics(a3(:,:,end), s3.ROI_info_table(1,1));
end

function [f1, f2, f3] = get_2nd_roi_values(a1, a2, a3, s1, s2, s3)

    if(numel(s1.ROI_info_table)==2)
        f1 = roi_statistics(a1(:,:,end), s1.ROI_info_table(1,1));
    else
        f1.m = -1;
    end
    
    if(numel(s2.ROI_info_table)==2)
        f2 = roi_statistics(a2(:,:,end), s2.ROI_info_table(1,1));
    else
        f2.m = -1;
    end
    
    if(numel(s3.ROI_info_table)==2)
        f3 = roi_statistics(a3(:,:,end), s3.ROI_info_table(1,1));
    else
        f3.m = -1;
    end
end
