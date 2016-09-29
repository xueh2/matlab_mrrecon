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
    
    copyfile(['D:\Temp\record_' GT_PORT '.txt'], [resDir '\record_' GT_HOST '_' GT_PORT '.txt']);
else
    [key, user] = sshKeyLookup(GT_HOST);
    if (~isempty(user) & startRemoteGT)
        StopGadgetronOnRemote(GT_HOST, GT_PORT);
        CopyGadgetronRecordOnRemote(GT_HOST, GT_PORT, [resDir '\record_' GT_HOST '_' GT_PORT '.txt']);
    end
end

cd(currDir)
