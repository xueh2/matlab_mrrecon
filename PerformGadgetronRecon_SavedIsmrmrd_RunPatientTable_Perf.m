
function tUsed = PerformGadgetronRecon_SavedIsmrmrd_RunPatientTable_Perf(dataDir, PerfTable, gt_host, resDir, cleanRemote, startRemoteGT, configName_preset, styleSheet)
% tUsed = PerformGadgetronRecon_SavedIsmrmrd_RunPatientTable_Perf(dataDir, PerfTable, gt_host, resDir, cleanRemote, startRemoteGT, configName_preset, styleSheet)
% tUsed = PerformGadgetronRecon_SavedIsmrmrd_RunPatientTable_Perf('I:\KAROLINSKA', PerfTable, 'localhost', 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area_Aif_recorded')
% tUsed = PerformGadgetronRecon_SavedIsmrmrd_RunPatientTable_Perf('I:\BARTS', PerfTable, 'localhost', 'D:\data\ut\NewData\PaperResults\BARTS_Area_Aif_recorded')
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

setenv('GT_HOST', gt_host); setenv('GT_PORT', GT_PORT);

if(nargin<4)
    resDir = dataDir;
end

if(nargin<5)
    cleanRemote = 0;
end

if(nargin<6)
    startRemoteGT = 1;
end

if(nargin<7)
    configName_preset = [];
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
if ( nargin >= 4 )
    styleSheetDefault = [ '%GADGETRON_DIR%\install\schema/' styleSheet];
end

xmlUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xml';

if(cleanRemote)
    [key, user] = sshKeyLookup(gt_host);
    gt_command = ['rm -rf /tmp/gadgetron_data/*'];
    command = ['ssh -i ' key ' ' user '@' gt_host ' "' gt_command '"'];
    command
    dos(command, '-echo');    
end

num = size(PerfTable, 1)-1;
tUsed = [];
for n=1:num
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, 18} ' - ' PerfTable{n+1, 19}]);   
    [tU, ig] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, PerfTable{n+1, 18}, gt_host, resDir, 0, 0, startRemoteGT, configName_preset, styleSheet);
    if(~isempty(tU))
        tUsed = [tUsed; {n, PerfTable{n+1, 18}, tU{2}, tU{3}}];
    end
    
    [tU, ig] = PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, PerfTable{n+1, 19}, gt_host, resDir, 0, 0, startRemoteGT, configName_preset, styleSheet);
    if(~isempty(tU))
        tUsed = [tUsed; {n, PerfTable{n+1, 19}, tU{2}, tU{3}}];
    end
end
