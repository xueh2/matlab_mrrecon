function t_linear = PerformGadgetronRecon_Perfusion_SavedIsmrmrd(dataName, gt_host, res_dir, xml, copy_debug, gt_port, startGT)
% t_linear = PerformGadgetronRecon_Perfusion_SavedIsmrmrd(dataName, gt_host, res_dir, xml, copy_debug, gt_port, startGT)

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

startRemoteGT = 1;

if(nargin>=6)
    GT_PORT = gt_port;
    GT_PORT
end

if(nargin>=7)
    startRemoteGT = startGT;
    startRemoteGT
end

setenv('GT_HOST', gt_host); setenv('GT_PORT', GT_PORT);

%% run the data
cd(res_dir)
command = ['gadgetron_ismrmrd_client -f ' dataName ' -C ' xml ' -a ' gt_host ' -p %GT_PORT% -F %OutputFormat% ']
tic; dos(command); t_linear = toc;

%% start process
UTDir = getenv('GTPLUS_UT_DIR')
home_ut = home(length(UTDir)+2:end)
UTCase_Flow = {home_ut, '20140408_13h17m13s_10066', 'VD',       xml,        linear_res_dir,              'grappa_flow_ref_PSIR',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL_Flow = {home_ut, '20140408_13h17m13s_10066', 'VD',    xmlNL,      nl_res_dir,              'slep_flow_ref_new',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL_Flow_Cloud = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping_NonLinear_Gateway.xml',        'slep_flow_cloud_res_new',              'slep_flow_cloud_ref_new',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}

% ------------------------------------------------------------

debugFolder = 'E:\gtuser\mrprogs\install\DebugOutput';

isVD = 1;
isVD11 = 1;

xslName = 'IsmrmrdParameterMap_Siemens_Perfusion.xsl';

%% for linear perf flow mapping
try
    rmdir(debugFolder, 's');
catch
end

home = fullfile(UTDir, home)
home_ut = home(length(UTDir)+2:end)
UTCase_Flow{1} = home_ut;

try
    mkdir(debugFolder);
catch
end

inds = strfind(subdir, '\');
v = subdir;
if(~isempty(inds))
    folderName = v(1:inds(end)-1)
else
    folderName = []
end
UTCase_Flow{1,1} = fullfile(home_ut, folderName)   

if(~isempty(inds))
    UTCase_Flow{1,2} = v(inds(end)+1:end)   
else
    UTCase_Flow{1,2} = v   
end

[names_dat, num_dat] = findFILE(fullfile(home, subdir), '*.dat')
h5Only = 0;

if(num_dat==0)
    [names_dat, num_dat] = findFILE(fullfile(home, subdir), '*.h5');
    h5Only = 1;
end

t_linear = 0;
t_nl = 0;

if(num_dat==0)
    return;
end

if(runLinear)
    fullPath = fullfile(home, subdir, UTCase_Flow{1, 5})
    command = ['rmdir /S /Q ' fullfile(fullPath, 'DebugOutput')];
    dos(command, '-echo');

    t_linear = run_gt_recon_case(UTCase_Flow{1, 2}, UTCase_Flow{1, 4}, UTCase_Flow, deleteh5, startRemoteGT, [], h5Only);    

    if(copy_debug)
        if(strcmp(getenv('GT_HOST'), 'localhost')==1)
            movefile(debugFolder, fullPath, 'f');
        else
            [key, user] = sshKeyLookup(getenv('GT_HOST'));
            debug_folder = ['/home/' user '/Debug/DebugOutput']
            CopyGadgetronDebugOutputOnRemote(getenv('GT_HOST'), debug_folder, fullPath, 1)
        end
    end
end

% ------------------------------

if(runNL)
    UTCase_NL_Flow{1,1} = fullfile(home_ut, folderName)   
    UTCase_NL_Flow{1, 2} = UTCase_Flow{1,2};

    fullPath = fullfile(home, subdir, UTCase_NL_Flow{1, 5})
    command = ['rmdir /S /Q ' fullfile(fullPath, 'DebugOutput')];
    dos(command, '-echo');

    t_nl = run_gt_recon_case(UTCase_NL_Flow{1, 2}, UTCase_NL_Flow{1, 4}, UTCase_NL_Flow, deleteh5, startRemoteGT, [], h5Only);    

    if(copy_debug)
        if(strcmp(getenv('GT_HOST'), 'localhost')==1)
            movefile(debugFolder, fullPath, 'f');
        else
            [key, user] = sshKeyLookup(getenv('GT_HOST'));
            debug_folder = ['/home/' user '/Debug/DebugOutput']
            CopyGadgetronDebugOutputOnRemote(getenv('GT_HOST'), debug_folder, fullPath, 1)
        end
    end
end

