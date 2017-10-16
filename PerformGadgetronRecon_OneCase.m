function t = PerformGadgetronRecon_OneCase(UT_Dir, home, subdir, gt_host, res_dir, deleteh5, xml, copy_debug, gt_port, startGT)
% t = PerformGadgetronRecon_OneCase(UT_Dir, home, subdir, gt_host, res_dir, deleteh5, xml, copy_debug, gt_port, startGT)
% PerformGadgetronRecon_OneCase('D', 'binning\CNMC_20160822', 'meas_MID00156_FID28880_AX_R2_RT_Cine_Binning_Res160_LocalCloud', 'barbados', 'res', deleteh5, 'GTPrep_2DT_RTCine_KspaceBinning_Cloud.xml', copy_debug, gt_port, startGT) 
% PerformGadgetronRecon_OneCase('D', 'binning\CNMC_20160822', 'meas_MID00041_FID28765_R2_RT_Cine_Binning_Res192_LocalCloud', 'barbados', 'res', deleteh5, 'Generic_Cartesian_NonLinear_Spirit_RealTimeCine_Cloud.xml', copy_debug, gt_port, startGT) 
% PerformGadgetronRecon_OneCase('T', 'binning\CNMC_20160822', 'meas_MID00041_FID28765_R2_RT_Cine_Binning_Res192_LocalCloud', 'barbados', 'res', deleteh5, 'Generic_Cartesian_NonLinear_Spirit_RealTimeCine_Cloud.xml', copy_debug, gt_port, startGT) 

if(nargin<5)
    res_dir = 'res';
end

if(nargin<6)
    deleteh5 = 0;
end

if(nargin<8)
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

if(nargin>=9)
    GT_PORT = gt_port;
    GT_PORT
end

if(nargin>=10)
    startRemoteGT = startGT;
    startRemoteGT
end

setenv('GT_HOST', gt_host); setenv('GT_PORT', GT_PORT);
getenv('GT_HOST')
getenv('GT_PORT')

%% start process
UTDir = getenv('GTPLUS_UT_DIR')
home_ut = home(length(UTDir)+2:end)
UTCase_Flow = {home_ut, '20140408_13h17m13s_10066', 'VD',       xml,        res_dir,              'grappa_flow_ref_PSIR',   'IsmrmrdParameterMap_Siemens.xsl'   ;}

% ------------------------------------------------------------

debugFolder = 'D:\gtuser\mrprogs\install\DebugOutput';

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

t = 0;

% if(num_dat==0)
%     return;
% end

fullPath = fullfile(home, subdir, UTCase_Flow{1, 5})
command = ['rmdir /S /Q ' fullfile(fullPath, 'DebugOutput')];
dos(command, '-echo');

t = run_gt_recon_case(UTCase_Flow{1, 2}, UTCase_Flow{1, 4}, UTCase_Flow, deleteh5, startRemoteGT, [], h5Only);    

if(copy_debug)
    if(strcmp(getenv('GT_HOST'), 'localhost')==1)
        movefile(debugFolder, fullPath, 'f');
    else
        [key, user] = sshKeyLookup(getenv('GT_HOST'));
        debug_folder = ['/home/' user '/Debug/DebugOutput']
        CopyGadgetronDebugOutputOnRemote(getenv('GT_HOST'), debug_folder, fullPath, 1)
    end
end

