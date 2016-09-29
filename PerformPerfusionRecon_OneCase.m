function [t_linear, t_nl] = PerformPerfusionRecon_OneCase(UT_Dir, home, subdir, gt_host, runLinear, runNL, linear_res_dir, nl_res_dir, deleteh5, xml, xmlNL, copy_debug, gt_port, startGT)
% [t_linear, t_nl] = PerformPerfusionRecon_OneCase(UT_Dir, home, subdir, gt_host, runLinear, runNL, linear_res_dir, nl_res_dir, deleteh5, xml, xmlNL, copy_debug, gt_port, startGT)
% PerformPerfusionRecon_OneCase('E', 'DualBolus\NewData\KAROLINSKA', subdir, 'palau', runLinear, runNL, 'grappa_flow_res_BTEX20_TwoCompWithShifts', 'slep_flow_res_BTEX20_TwoCompWithShifts', 0)
% PerformPerfusionRecon_OneCase('F', 'DualBolus\FFR', subdir, 'palau', runLinear, runNL, 'grappa_flow_res_BTEX20_TwoCompWithShifts', 'slep_flow_res_BTEX20_TwoCompWithShifts', 0)

if(nargin<7)
    linear_res_dir = 'grappa_flow_res_BTEX20_TwoCompWithShifts';
end

if(nargin<8)
    linear_res_dir = 'slep_flow_res_BTEX20_TwoCompWithShifts';
end

if(nargin<9)
    deleteh5 = 0;
end

if(nargin<10)
    xml = 'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping.xml';
end

if(nargin<11)
    xmlNL = 'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping_NonLinear.xml';
end

if(nargin<12)
    copy_debug = 0;
end

set_UT_Dir(UT_Dir)

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

if(nargin>=13)
    GT_PORT = gt_port;
    GT_PORT
end

if(nargin>=14)
    startRemoteGT = startGT;
    startRemoteGT
end

setenv('GT_HOST', gt_host); setenv('GT_PORT', GT_PORT);
getenv('GT_HOST')
getenv('GT_PORT')

%% start process
UTDir = getenv('GTPLUS_UT_DIR')
home_ut = home(length(UTDir)+2:end)
UTCase_Flow = {home_ut, '20140408_13h17m13s_10066', 'VD',       xml,        linear_res_dir,              'grappa_flow_ref_PSIR',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL_Flow = {home_ut, '20140408_13h17m13s_10066', 'VD',    xmlNL,      nl_res_dir,              'slep_flow_ref_new',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL_Flow_Cloud = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping_NonLinear_Gateway.xml',        'slep_flow_cloud_res_new',              'slep_flow_cloud_ref_new',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}

% ------------------------------------------------------------

debugFolder = 'D:\gtuser\mrprogs\install\DebugOutput';

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

