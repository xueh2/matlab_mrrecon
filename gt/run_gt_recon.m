function timeUsed = run_gt_recon(folderDir, dataName, h5Name, deleteh5, isVD, isVD11, isAdjScan, configName, resDir, styleSheet, startRemoteGT, h5Only, remoteXml, compressionBit)
% run_gt_recon(folderDir, dataName, h5Name, deleteh5, isVD, isVD11, isAdjScan, configName, resDir, styleSheet, startRemoteGT, h5Only)

% siemens to hdf5

currDir = pwd;

GTHome = getenv('GADGETRON_HOME')
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
if(remoteXml)
    configNameUsed = configName;
else
    configNameUsed = fullfile(GTConfigFolder, configName);
end

GT_HOST = getenv('GT_HOST')
GT_PORT = getenv('GT_PORT')

if(startRemoteGT)
    if((strcmp(GT_HOST, 'localhost')==1))
        cd('D:\gtuser\gt_windows_setup')
        command = ['gadgetron -p %GT_PORT% > D:\Temp\record_' GT_PORT '.txt']
        dos([command ' &'])
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

styleSheetUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl';
if ( nargin >= 8 )
    styleSheetUsed = [ '%GADGETRON_DIR%\install\schema/' styleSheet];
end

flashRef = 0;
if( ~isempty(strfind(styleSheet, 'IsmrmrdParameterMap_Siemens_EPI_FLASHREF')) )
    flashRef = 1;
end

xmlUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xml';

% if VD, run the dependent scan

lenC = length(configName);
currPath = pwd;

lenUsed = 200 - (1+length(currPath));   
if(lenUsed>lenC)
    lenUsed = lenC;
end

configNameShortened = configName(1:lenUsed);

date_suffix = datestr(date, 'yyyymmdd');

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

if ( isVD )
    if ( isVD11 )
        if ( ~isFileExist(h5Name) || deleteh5 )
            delete(h5Name);
            command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 -X -H --studyDate ' studyDate]
            dos(command);
            
            command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 --studyDate ' studyDate]
            dos(command);
        end
        
        command = ['gadgetron_ismrmrd_client -f ' h5Name xml_opt configNameUsed ' -a %GT_HOST% -p %GT_PORT% -F %OutputFormat% -G ' configNameShortened ' -o ref_' date_suffix '.h5' compression_opt]
        command
        tic; dos(command); timeUsed = toc;
    else

        if(h5Only)
            
            [names, num] = findFILE(folderDir, 'ISMRMRD_Noise_dependency_*.h5');          
            for n=1:num
                command = ['gadgetron_ismrmrd_client -f ' names{n} ' -c default_measurement_dependencies.xml -a %GT_HOST% -p %GT_PORT% ']
                dos(command, '-echo');
            end            
            
            command = ['gadgetron_ismrmrd_client -f ' h5Name xml_opt configNameUsed ' -a %GT_HOST% -p %GT_PORT% -F %OutputFormat% -G ' configNameShortened ' -o ref_' date_suffix '.h5' compression_opt]
            command
            tic; dos(command); timeUsed = toc;
        else
            hasAdj = 1;
            if(deleteh5)
                % check whether two measurements are in the dat file
                command = ['siemens_to_ismrmrd -f ' dataName ' -o test_meas.h5 --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 --studyDate ' studyDate]
                dos(command, '-echo');

                if(~isFileExist('test_meas.h5'))
                    hasAdj = 0;
                else
                    delete('test_meas.h5');
                end
            end

            [names, num] = findFILE(folderDir, '*AdjCoilSens*.dat');          
            for n=1:num
                command = ['siemens_to_ismrmrd -f ' names{n} ' -o ' names{n} '_AdjCoilSens.h5 --user-map %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xml --user-stylesheet %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl -z 1 --studyDate ' studyDate]
                dos(command, '-echo');
                
                command = ['gadgetron_ismrmrd_client -f ' names{n} '_AdjCoilSens.h5 -c default_measurement_dependencies.xml -a %GT_HOST% -p %GT_PORT% ']
                dos(command, '-echo');
            end            
                       
            if(hasAdj)
                if ( ~isFileExist('AdjCoilSens.h5') || deleteh5 )
                    delete('AdjCoilSens.h5');
                    command = ['siemens_to_ismrmrd -f ' dataName ' -o AdjCoilSens.h5 --user-map %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xml --user-stylesheet %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl -z 1 --studyDate ' studyDate]
                    dos(command, '-echo');
                end

                command = ['gadgetron_ismrmrd_client -f AdjCoilSens.h5 -c default_measurement_dependencies.xml -a %GT_HOST% -p %GT_PORT% ']
                dos(command, '-echo');
            end

            if ( ~isFileExist(h5Name) || deleteh5 )

                delete(h5Name);

                if(hasAdj)                            
                    command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 -X -H --studyDate ' studyDate]
                    if(flashRef)
                        command = [command ' -F'];
                    end
                    dos(command, '-echo');

                    command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 --studyDate ' studyDate]
                    if(flashRef)
                        command = [command ' -F'];
                    end
                    dos(command, '-echo');
                else
                    command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 1 -X -H --studyDate ' studyDate]
                    if(flashRef)
                        command = [command ' -F'];
                    end
                    dos(command, '-echo');

                    command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 1 --studyDate ' studyDate]
                    if(flashRef)
                        command = [command ' -F'];
                    end
                    dos(command, '-echo');
                end
            end

            command = ['gadgetron_ismrmrd_client -f ' h5Name xml_opt configNameUsed ' -a %GT_HOST% -p %GT_PORT% -F %OutputFormat% -G ' configNameShortened ' -o ref_' date_suffix '.h5' compression_opt]
            command
            tic; dos(command); timeUsed = toc;
        end
    end
else  
    if(isAdjScan)
        command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xml --user-stylesheet %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl -z 1 -X -H --studyDate ' studyDate]
        dos(command);
            
        command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xml --user-stylesheet %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl -z 1 -X --studyDate ' studyDate]
        dos(command);
            
        command = ['gadgetron_ismrmrd_client -f ' h5Name ' -c default_measurement_dependencies.xml -a %GT_HOST% -p %GT_PORT% ']
        tic; dos(command); timeUsed = toc;
    else
        if ( ~isFileExist(h5Name) || deleteh5 )
            delete(h5Name);
            command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_VB17.xml --user-stylesheet ' styleSheetUsed ' -z 1 -X -H --studyDate ' studyDate]
            dos(command);

            command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_VB17.xml --user-stylesheet ' styleSheetUsed ' -z 1 --studyDate ' studyDate]
            dos(command);
        end

        command = ['gadgetron_ismrmrd_client -f ' h5Name xml_opt configNameUsed ' -a %GT_HOST% -p %GT_PORT% -F %OutputFormat% -G ' configNameShortened ' -o ref_' date_suffix '.h5' compression_opt]
        command
        tic; dos(command); timeUsed = toc;
    end
end

if(strcmp(GT_HOST, 'localhost')==1)
    if(startRemoteGT)
        command = ['taskkill /F /FI "IMAGENAME eq gadgetron.*"'];
        dos(command)
    end
    
    try
        copyfile(['D:\Temp\record_' GT_PORT '.txt'], [resDir '\record_' GT_HOST '_' GT_PORT '.txt']);
    catch
    end
else
    [key, user] = sshKeyLookup(GT_HOST);
    if (~isempty(user) & startRemoteGT)
        StopGadgetronOnRemote(GT_HOST, GT_PORT);
        CopyGadgetronRecordOnRemote(GT_HOST, GT_PORT, [resDir '\record_' GT_HOST '_' GT_PORT '.txt']);
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
