
function [tUsed, ignored] = PerformGadgetronRecon_WholeStudy(dataDir, gt_host, startRemoteGT, deleteh5, start_case, end_case, ignore_list, styleSheet)
% perform gadgetron reconstruction for the whole study
% tUsed = PerformGadgetronRecon_WholeStudy(dataDir, gt_host, startRemoteGT, deleteh5, start_case, end_case, ignore_list, styleSheet)

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

if(nargin<3)
    startRemoteGT = 0;
end

if(nargin<4)
    deleteh5 = 0;
end

if(nargin<5)
    start_case = 1;
end

if (nargin<6)
    end_case = 1000;
end

if(nargin<7)
    ignore_list = {'retrogated', 'mini', 'test_', 'tag', 'qs', 'qp'};
end

if(nargin<8)
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

% ------------------------------------------------------------

% find all dependent scan and process them first
[names, num] = findFILE(dataDir, '*AdjCoilSens*.dat');          
for n=1:num
    
    [pathstr, name, ext] = fileparts(names{n});
    
    mkdir(fullfile(dataDir, name));
    cd(fullfile(dataDir, name))
    
    h5Name = [names{n} '_AdjCoilSens.h5'];
    
    if(~isFileExist(h5Name) | deleteh5)
        delete(h5Name);
        command = ['siemens_to_ismrmrd -f ' names{n} ' -o ' h5Name ' --user-map %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xml --user-stylesheet %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl -z 1 ']
        dos(command, '-echo');
    end
    
    command = ['gadgetron_ismrmrd_client -f ' h5Name ' -c default_measurement_dependencies.xml -a %GT_HOST% -p %GT_PORT% ']
    dos(command, '-echo');
end 

tUsed = [];
ignored = [];

[names, num_dat] = findFILE(dataDir, '201*.dat'); 

auto_saved = 1;

if(num==0 & num_dat>0)
    disp('Auto saved study detected ...');
    [names, num] = findFILE(dataDir, '201*.dat');
    auto_saved = 1;
    hasAdj = 1;
else
    disp('Manually saved study detected ...');   
    auto_saved = 0;    
    hasAdj = 0;    
    [names, num] = findFILE(dataDir, 'meas_*.dat');               
end

for n=1:num
    disp([num2str(n) ' - ' names{n}])
end

% run the recon
if(end_case>num)
    end_case = num;
end

for n=start_case:end_case

    [pathstr, name, ext] = fileparts(names{n});

    if(~auto_saved & ~isempty(strfind(name, 'Adj')))
        continue;
    end

    mkdir(fullfile(dataDir, name));
    cd(fullfile(dataDir, name))

    delete(fullfile(dataDir, name, 'res*.h5'));
    delete(fullfile(dataDir, name, 'out*.h5'));
    delete(fullfile(dataDir, name, '*.xml'));

    delete(fullfile(dataDir, name, '*.nii'));
    delete(fullfile(dataDir, name, 'gadgetron_*.hdr'));
    delete(fullfile(dataDir, name, 'gadgetron_*.img'));
    delete(fullfile(dataDir, name, 'Generic*.hdr'));
    delete(fullfile(dataDir, name, 'Generic*.img'));
    delete(fullfile(dataDir, name, 'GTPrep*.hdr'));
    delete(fullfile(dataDir, name, 'GTPrep*.img'));
    delete(fullfile(dataDir, name, 'GT*.hdr'));
    delete(fullfile(dataDir, name, 'GT*.img'));
    delete(fullfile(dataDir, name, '*.attrib'));

    if(deleteh5)
        delete(fullfile(dataDir, name, '*.xml'));
        delete(fullfile(dataDir, name, '*.h5'));
    end
        
    dataName = names{n};
    h5Name = fullfile(dataDir, name, [name '.h5']);

    if(deleteh5)
        delete(h5Name);
    end

    convertsionFlag = [];
        
    % handle dependencies
    if(auto_saved)
        if(~isFileExist(h5Name))      
            command = ['siemens_to_ismrmrd -f ' names{n} ' -o ' names{n} '_AdjCoilSens.h5 --user-map %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xml --user-stylesheet %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl -z 1 ']
            dos(command, '-echo');

            command = ['gadgetron_ismrmrd_client -f ' names{n} '_AdjCoilSens.h5 -c default_measurement_dependencies.xml -a %GT_HOST% -p %GT_PORT% ']
            dos(command, '-echo');
        end 
    else
        if(~isFileExist(h5Name))               
            if(hasAdj)
                command = ['siemens_to_ismrmrd -f ' dataName ' -o ' names{n} '_AdjCoilSens.h5 --user-map %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xml --user-stylesheet %GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens.xsl -z 1 ']
                dos(command, '-echo');

                command = ['gadgetron_ismrmrd_client -f ' dataName '_AdjCoilSens.h5 -c default_measurement_dependencies.xml -a %GT_HOST% -p %GT_PORT% ']
                dos(command, '-echo');
            end
        end
    end
    
    % get the data prot
    styleSheetUsed = styleSheetDefault;
        
    if(auto_saved)
        if(~isFileExist(h5Name))      
            command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 -X -H ' convertsionFlag]
            dos(command);
            delete(h5Name)
        end 
    else
        if(~isFileExist(h5Name))               
            if(hasAdj)
                command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 -X -H ' convertsionFlag]
                dos(command);
                delete(h5Name)
            else
                command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 1 -X -H ' convertsionFlag]
                dos(command);
                delete(h5Name)
            end
        end
    end
    
    % find the Ice ipr name
    configName = findGadgetronXmlConfigName(fullfile(dataDir, name), dataName)
    if(~isempty(strfind(configName, '_EPI')))
        convertsionFlag = '-F';
    end
    
    %% ignore some scans
    ignore_flag = 0;
    for kk=1:numel(ignore_list)
        if(~isempty(strfind(lower(configName), ignore_list{kk})) | ~isempty(strfind(lower(name), ignore_list{kk})))
            ignore_flag = 1;
            break;
        end
    end
    
    finfo = dir(dataName);
    
    if(ignore_flag)
        disp(['Ignore - ' num2str(n) ' - ' name]);
        
        cd(dataDir)
        rmdir(fullfile(dataDir, name), 's');
        
        ignored = [ignored; {n, name, finfo.bytes/1024}];
        
        continue;
    end
        
    if(finfo.bytes<2*2048*2048)
        disp(['File size too small - ' num2str(n) ' - ' name]);
        cd(dataDir)
        delete(fullfile(dataDir, name, '*.*'));
        rmdir(fullfile(dataDir, name), 's');
        ignored = [ignored; {n, name, finfo.bytes/1024}];
        continue;
    end
    
    %% run the scan
    perfXSL = 0;
    if(~isempty(strfind(configName, 'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2')))
        styleSheetUsed = styleSheetPerfusionUsed;
        perfXSL = 1;
    end
    
    if(auto_saved)
        if(~isFileExist(h5Name))             
            if(perfXSL)
                command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 -X -H ' convertsionFlag]
                dos(command);
            end
            
            command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 ' convertsionFlag]
            dos(command);
        end 
    else
        if(~isFileExist(h5Name))               
            if(hasAdj)              
                if(perfXSL)
                    command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 -X -H ' convertsionFlag]
                    dos(command);
                end
                
                command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 2 ' convertsionFlag]
                dos(command);
            else
                if(perfXSL)
                    command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 1 -X -H ' convertsionFlag]
                    dos(command);
                end
                
                command = ['siemens_to_ismrmrd -f ' dataName ' -o ' h5Name ' --user-map ' xmlUsed ' --user-stylesheet ' styleSheetUsed ' -z 1 ' convertsionFlag]
                dos(command);
            end
        end
    end
    
    convertsionFlag = [];
    
    if (~isempty(configName))
        % perform recon
        configNameUsed = fullfile(GTConfigFolder, configName);

        lenC = length(configName);
        currPath = pwd;

        lenUsed = 200 - (1+length(currPath));   
        if(lenUsed>lenC)
            lenUsed = lenC;
        end

        configNameShortened = configName(1:lenUsed);

        command = ['gadgetron_ismrmrd_client -f ' h5Name ' -C ' configNameUsed ' -a %GT_HOST% -p %GT_PORT% -F %OutputFormat% -G ' configNameShortened ' -o ref_' date_suffix '.h5']
        tic; dos(command); timeUsed = toc;
            
        tUsed = [tUsed; {n, name, timeUsed, configName}];
    end
end
