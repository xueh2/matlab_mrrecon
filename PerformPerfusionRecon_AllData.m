
clear all
close all

% ------------------------------------------------------------

if(0)
    setenv('GT_HOST', 'denmark'); setenv('GT_PORT', '9006');
    setenv('GT_HOST', 'barbados'); setenv('GT_PORT', '9006');
    setenv('GT_HOST', 'palau'); setenv('GT_PORT', '9008');
    setenv('GT_HOST', 'localhost'); setenv('GT_PORT', '9002');
    setenv('GT_HOST', 'samoa'); setenv('GT_PORT', '9016');
    startRemoteGT = 1;
end

getenv('GT_HOST')
getenv('GT_PORT')

% ------------------------------------------------------------
home_all = [];

set_UT_Dir('E')
UTDir = getenv('GTPLUS_UT_DIR')
homes = {'DualBolus\FLASH', 'DualBolus\SSFP', 'DualBolus\SUB_Perf\NewSeq', 'DualBolus\FFR', ... 
    'DualBolus\NewData\FLASH', 'DualBolus\NewData\SSFP', 'DualBolus\NewData\SUB_1p5T', 'DualBolus\NewData\KAROLINSKA'};

% setenv('GTPLUS_UT_DIR', 'F:\data');
% UTDir = getenv('GTPLUS_UT_DIR')
% homes = {'DualBolus\highSpatial', 'DualBolus\highTemporal', 'DualBolus\slowPerf', 'DualBolus\TwoInjectionRate'};

for ii=1:numel(homes)
    home = fullfile(UTDir, homes{ii});
    home_all = [home_all; {home}];
end
home_all

%% start process

home_ut = home(length(UTDir)+2:end)
UTCase = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2.xml', 'grappa_res',              'grappa_ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_NonLinear.xml',           'slep_res',          'slep_ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL_Cloud = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_NonLinear_Gateway.xml',           'slep_cloud_res',          'slep_cloud_ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}

UTCase_Flow = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping.xml',                                   'grappa_flow_res_BTEX20_TwoCompWithShifts',              'grappa_flow_ref_PSIR',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL_Flow = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping_NonLinear.xml',                      'slep_flow_res_BTEX20_TwoCompWithShifts',              'slep_flow_ref_new',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL_Flow_Cloud = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping_NonLinear_Gateway.xml',        'slep_flow_cloud_res_new',              'slep_flow_cloud_ref_new',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}

% ------------------------------------------------------------

exDir = {'ICERecon', 'grappa', 'TXMapping', 'moco', 'mocoSyn', 'mocoPS' , 'mocoPSSyn', 'seg', 'magFitting', 'slep_res', 'slep_cloud_res', 'grappa_res', 'ICE', 'DebugOutput', 'PhantomT1maps', 'slep_flow_res', 'grappa_flow_res', 'slep_cloud_flow_res', 'slep_cloud_flow_res2'};
[subdir, num] = FindAllEndDirectoryExclusive(home, exDir)

debugFolder = 'E:\gtuser\mrprogs\install\DebugOutput';

isVD = 1;
isVD11 = 1;

configName = 'GT_2DT_Perfusion_NonLinear.xml';
res = 'grappa_res';

xslName = 'IsmrmrdParameterMap_Siemens_Perfusion.xsl';
deleteh5 = 0;
% deleteh5 = 1;

runLinear = 1;
runNL = 1;
% runNL = 0;

%% for linear perf flow mapping
rmdir(debugFolder, 's');

% for pp=4:numel(home_all)
for pp=1:numel(home_all)
    
    home = home_all{pp}
    
%     [subdir, num] = FindAllEndDirectoryExclusive(home, exDir)
    
    home_ut = home(length(UTDir)+2:end)
    UTCase_Flow{1} = home_ut;

    % ------------------------------------------------------------

    exDir = {'ICERecon', 'grappa', 'TXMapping', 'moco', 'mocoSyn', 'mocoPS' , 'mocoPSSyn', 'seg', 'magFitting', 'slep_res', 'slep_cloud_res', 'grappa_res', 'ICE', 'DebugOutput', 'PhantomT1maps', 'slep_flow_res', 'grappa_flow_res', 'slep_cloud_flow_res', 'slep_cloud_flow_res2', 'T2Map'};
    [subdir, num] = FindAllEndDirectoryExclusive(home, exDir)

    % ------------------------------------------------------------
    
    startI = 1;
    
    t_linear = [];
    t_nl = [];
    
    for ii=startI:num  

        if ( ~isempty(strfind( lower(subdir{ii}), 'mini')) | ~isempty(strfind(subdir{ii}, 'MINI')) | ~isempty(strfind(subdir{ii}, 'Mini')) )
            continue;
        end
        
        if ( ~isempty(strfind( lower(subdir{ii}), 't1map')) )
            continue;
        end
       
%         if ( isempty(strfind( lower(subdir{ii}), 'rest')) &  isempty(strfind( lower(subdir{ii}), 'stress')) )
%             continue;
%         end
        
        mkdir(debugFolder);
   
        inds = strfind(subdir{ii}, '\');
        v = subdir{ii};
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

        [names_dat, num_dat] = findFILE(fullfile(home, subdir{ii}), '*.dat')
        h5Only = 0;
        
        if(num_dat==0)
            [names_dat, num_dat] = findFILE(fullfile(home, subdir{ii}), '*.h5');
            h5Only = 1;
        end
        
        if(num_dat==0)
            continue;
        end
        
        if(runLinear)
            fullPath = fullfile(home, subdir{ii}, UTCase_Flow{1, 5})
            command = ['rmdir /S /Q ' fullfile(fullPath, 'DebugOutput')];
            dos(command, '-echo');

            t = run_gt_recon_case(UTCase_Flow{1, 2}, UTCase_Flow{1, 4}, UTCase_Flow, deleteh5, startRemoteGT, h5Only);    

            t_linear = [t_linear; {subdir{ii} t}];
            
            if(strcmp(getenv('GT_HOST'), 'localhost')==1)
                movefile(debugFolder, fullPath, 'f');
            else
                [key, user] = sshKeyLookup(getenv('GT_HOST'));
                debug_folder = ['/home/' user '/Debug/DebugOutput']
                CopyGadgetronDebugOutputOnRemote(getenv('GT_HOST'), debug_folder, fullPath, 1)
            end 
        end
        
        % ------------------------------
        
        if(runNL)
            UTCase_NL_Flow{1,1} = fullfile(home_ut, folderName)   
            UTCase_NL_Flow{1, 2} = UTCase_Flow{1,2};

            fullPath = fullfile(home, subdir{ii}, UTCase_NL_Flow{1, 5})
            command = ['rmdir /S /Q ' fullfile(fullPath, 'DebugOutput')];
            dos(command, '-echo');

            t = run_gt_recon_case(UTCase_NL_Flow{1, 2}, UTCase_NL_Flow{1, 4}, UTCase_NL_Flow, deleteh5, startRemoteGT, h5Only);    

            t_nl = [t_nl; {subdir{ii} t}];
            
            if(strcmp(getenv('GT_HOST'), 'localhost')==1)
                movefile(debugFolder, fullPath, 'f');
            else
                [key, user] = sshKeyLookup(getenv('GT_HOST'));
                debug_folder = ['/home/' user '/Debug/DebugOutput']
                CopyGadgetronDebugOutputOnRemote(getenv('GT_HOST'), debug_folder, fullPath, 1)
            end 
        end
    end
    
%     for ii=1:num  
%         inds = strfind(subdir{ii}, '\');
%         v = subdir{ii};
%         if(~isempty(inds))
%             folderName = v(1:inds(end)-1)
%         else
%             folderName = []
%         end
%         src = fullfile(home_ut, folderName) 
%         
%         if(~isempty(inds))
%             data = v(inds(end)+1:end)   
%         else
%             data = v   
%         end
%         
%         dstDir = fullfile(dstRoot, src, data)
%         mkdir(dstDir);
%         copyfile(fullfile(home, folderName, data, [data '.dat']), dstDir);        
%     end

%     for ii=startI:num  
% 
%         if ( ~isempty(strfind( lower(subdir{ii}), 'mini')) | ~isempty(strfind(subdir{ii}, 'MINI')) | ~isempty(strfind(subdir{ii}, 'Mini')) )
%             continue;
%         end
%         
%         if ( ~isempty(strfind( lower(subdir{ii}), 't1map')) )
%             continue;
%         end
%                 
%         inds = strfind(subdir{ii}, '\');
%         v = subdir{ii};
%         if(~isempty(inds))
%             folderName = v(1:inds(end)-1)
%         else
%             folderName = []
%         end
%         UTCase_NL_Flow{1,1} = fullfile(home_ut, folderName)   
% 
%         if(~isempty(inds))
%             UTCase_NL_Flow{1,2} = v(inds(end)+1:end)   
%         else
%             UTCase_NL_Flow{1,2} = v   
%         end
% 
%         [names_dat, num_dat] = findFILE(fullfile(home, subdir{ii}), '*.dat')
%         
%         if(num_dat==0)
%             continue;
%         end
%         
%         delete(fullfile(home, subdir{ii}, UTCase_NL_Flow{1, 5}, 'DebugOutput/data_*'));
%         delete(fullfile(home, subdir{ii}, UTCase_NL_Flow{1, 5}, 'DebugOutput/workOrder_ref*'));
%         delete(fullfile(home, subdir{ii}, UTCase_NL_Flow{1, 5}, 'DebugOutput/ref_*'));
%     end

end

%% copy data
dstRoot = '\\137.187.134.135\share\DATA';

for pp=1:numel(home_all)
    
    home = home_all{pp}
    
    [subdir, num] = FindAllEndDirectoryExclusive(home, exDir)
    
    home_ut = home(length(UTDir)+2:end)
    UTCase_Flow{1} = home_ut;

    % ------------------------------------------------------------

    exDir = {'ICERecon', 'grappa', 'TXMapping', 'moco', 'mocoSyn', 'mocoPS' , 'mocoPSSyn', 'seg', 'magFitting', 'slep_res', 'slep_cloud_res', 'grappa_res', 'ICE', 'DebugOutput', 'PhantomT1maps', 'slep_flow_res', 'grappa_flow_res', 'slep_cloud_flow_res', 'slep_cloud_flow_res2'};
    [subdir, num] = FindAllEndDirectoryExclusive(home, exDir)

    for ii=1:num  
        inds = strfind(subdir{ii}, '\');
        v = subdir{ii};
        if(~isempty(inds))
            folderName = v(1:inds(end)-1)
        else
            folderName = []
        end
        src = fullfile(home_ut, folderName) 
        
        if(~isempty(inds))
            data = v(inds(end)+1:end)   
        else
            data = v   
        end
        
        dstDir = fullfile(dstRoot, src, data)
        mkdir(dstDir);
        copyfile(fullfile(home, folderName, data, [data '.dat']), dstDir);        
    end
end

%% for nonlinear recon
for ii=1:num  
       
    if ( ~isempty(strfind(subdir{ii}, 'mini')) | ~isempty(strfind(subdir{ii}, 'MINI'))  )
        continue;
    end
    
    inds = strfind(subdir{ii}, '\');
    v = subdir{ii};   
    if(~isempty(inds))
        folderName = v(1:inds(end)-1)
    else
        folderName = []
    end
    UTCase_NL_Flow{1,1} = fullfile(home_ut, folderName)   
    
    if(~isempty(inds))
        UTCase_NL_Flow{1,2} = v(inds(end)+1:end)   
    else
        UTCase_NL_Flow{1,2} = v   
    end
    
    t = run_gt_recon_case(UTCase_NL_Flow{1, 2}, UTCase_NL_Flow{1, 4}, UTCase_NL_Flow, deleteh5);    
   
end

%% for nonlinear cloud recon
for ii=1:num  
       
    inds = strfind(subdir{ii}, '\');
    v = subdir{ii};
    
    if(~isempty(inds))
        folderName = v(1:inds(end)-1)
    else
        folderName = []
    end
    UTCase_NL_Flow_Cloud{1,1} = fullfile(home_ut, folderName)   
    
    if(~isempty(inds))
        UTCase_NL_Flow_Cloud{1,2} = v(inds(end)+1:end)   
    else
        UTCase_NL_Flow_Cloud{1,2} = v   
    end
    
    t = run_gt_recon_case(UTCase_NL_Flow_Cloud{1, 2}, UTCase_NL_Flow_Cloud{1, 4}, UTCase_NL_Flow_Cloud, deleteh5);    
   
end

%% for non-linear perf flow mapping
rmdir(debugFolder, 's');
for ii=1:num 
      
    if ( ~isempty(strfind(subdir{ii}, 'mini')) | ~isempty(strfind(subdir{ii}, 'MINI'))  )
%         continue;
    end
    
    mkdir(debugFolder);
    mkdir(fullfile(debugFolder, '0'));
    mkdir(fullfile(debugFolder, '1'));
    mkdir(fullfile(debugFolder, '2'));
    mkdir(fullfile(debugFolder, '3'));
    mkdir(fullfile(debugFolder, '4'));
    mkdir(fullfile(debugFolder, '5'));
    mkdir(fullfile(debugFolder, '6'));
    mkdir(fullfile(debugFolder, '7'));
    
    inds = strfind(subdir{ii}, '\');
    v = subdir{ii};
    if(~isempty(inds))
        folderName = v(1:inds(end)-1)
    else
        folderName = []
    end
    UTCase_NL_Flow{1,1} = fullfile(home_ut, folderName)   
    
    if(~isempty(inds))
        UTCase_NL_Flow{1,2} = v(inds(end)+1:end)   
    else
        UTCase_NL_Flow{1,2} = v   
    end
   
    fullPath = fullfile(home, subdir{ii}, UTCase_NL_Flow{1, 5})
    commnd = ['rmdir /S /Q ' fullfile(fullPath, 'DebugOutput')];
    
    t = run_gt_recon_case(UTCase_NL_Flow{1, 2}, UTCase_NL_Flow{1, 4}, UTCase_NL_Flow, deleteh5);    
        
    movefile(debugFolder, fullPath, 'f');
end

%% image reviewing
for ii=1:num  
   
    subdir{ii}  

    cd(fullfile(home, subdir{ii}))
  
    try
        [flowmaps, acq_time, physio_time] = readGTPlusExportImageSeries('./grappa_flow_res', 115, 1);
        flowmaps = squeeze(flowmaps/1000);
        h = figure; imagescn(flowmaps(:,:,:,end), [0 6]); PerfColorMap;
        saveas(h, 'flowmaps', 'fig');
    catch
    end
    
    % if ( ~isFileExist(fullfile(home, subdir{ii}, 'DebugOutput', 'DebugOutput', 'aif_LV_mask.hdr')) )
    if ( ~isFileExist(fullfile(home, subdir{ii}, 'DebugOutput', 'aif_LV_mask.hdr')) )
        continue;
    end
    
    cd(fullfile(home, subdir{ii}, 'DebugOutput'))
    % cd(fullfile(home, subdir{ii}, 'DebugOutput', 'DebugOutput'))
       
    aif_moco_upsampled = analyze75read('aif_moco_upsampled.hdr');
    aif_moco_second_echo_upsampled = analyze75read('aif_moco_second_echo_upsampled.hdr');
    
    aif_mask = analyze75read('aif_LV_mask.hdr');
    aif_cin_all_echo0_signal = analyze75read('aif_cin_all_echo0_signal.hdr'); 
    aif_cin_all_echo0_signal_after_R2StarCorrection = analyze75read('aif_cin_all_echo0_signal_after_R2StarCorrection.hdr'); 
    aif_cin_all_echo1_signal = analyze75read('aif_cin_all_echo1_signal.hdr'); 
    
    aif_cin = analyze75read('aif_cin.hdr');
    Input_cin_computeFlowMap = analyze75read('Input_cin_computeFlowMap.hdr');
   
    aif_cin_all_echo0_LUTCorrection = analyze75read('aif_cin_echo0_all_signal_baseline_corrected.hdr');
    
    aif_cin_LUT_flash_pd_flash_sr = analyze75read('aif_cin_LUT_flash_pd_flash_sr.hdr');

    Input_perf_computeFlowMap_0 = analyze75read('Input_perf_computeFlowMap_0.hdr');
       
    cd(fullfile(home, subdir{ii}))
    
    mask_file = 'sr_mask.mat';
    if(~isFileExist(mask_file))
        figure; imagescn(Input_perf_computeFlowMap_0, [], [] ,[], 3);
        pause;
    end

    CASignal = roi_timeseries(Input_perf_computeFlowMap_0, mask_file, 1, 1);
    
    h = figure; imagescn(aif_mask);
    saveas(h, 'aif_mask', 'fig');
         
    h = figure; hold on; plot(aif_cin_all_echo0_signal); plot(aif_cin_all_echo0_signal_after_R2StarCorrection, 'b-.'); plot(aif_cin_all_echo1_signal, 'r'); hold off
    legend('aif echo 0', 'aif echo 1', 'aif echo 0 after R2* correctoin');
    title('aif R2* correction');
    saveas(h, 'aif_R2Star_correction', 'fig');
    
    h = figure; 
    hold on
    plot(aif_cin_all_echo0_LUTCorrection);

    v = mean(CASignal.m(1:5));
    plot(CASignal.m-v, 'r');
    hold off
    title('aif and SR in Gd');
    legend('aif', 'myo');
    saveas(h, 'aif_myo_in_Gd', 'fig');
    
    [path, name, ext] = fileparts(fullfile(home, subdir{ii}));
    
    save([name '_aif_record_new'], 'CASignal', 'aif_moco_upsampled', 'aif_moco_second_echo_upsampled', 'aif_mask', 'aif_cin_all_echo0_signal', 'aif_cin_all_echo1_signal', 'aif_cin_all_echo0_signal_after_R2StarCorrection', 'aif_cin_all_echo0_LUTCorrection', 'CASignal_LUTCorrection_0', 'aif_cin_LUT_flash_pd_flash_sr');
    
    closeall
end

%% make movie

series_ori = 103;
series = 104;
series_norm = 105;

folder = './grappa_res';
folder = './grappa_flow_res';


folder = './slep_flow_res';
series_ori = 103;
series = 104;
series_norm = 111;

PD = 3;

for ii =1:num
    
    dataDir = fullfile(home, subdir{ii}) 
    cd(dataDir)
    
%     ori = readGTPlusExportImageSeries_Squeeze(folder, series_ori);
%     
%     KLT = readGTPlusExportImageSeries_Squeeze(folder, series);
    
    KLT_norm = readGTPlusExportImageSeries_Squeeze(folder, series_norm);
    
    % ee = cat(5, ori(:,:,:,PD+1:end), KLT(:,:,:,PD+1:end));
    ee = KLT_norm;
    
    ee = flipdim(ee, 1);
    ee = flipdim(ee, 2);
    
    figure; imagescn(ee, [0 3000], [1 3], 8, 4);
    
    pause
    closeall
end


%% plot the myo of mini and aif of big

clear all
close all

home = 'D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150601\stress';
home = 'D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150518\stress';
home = 'D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150506\stress';

[subdir, num] = FindSubDirs(home)

mini = subdir{1}
full = subdir{2}

cd(fullfile(home, mini))

mini_im = readGTPlusExportImageSeries_Squeeze('./grappa_flow_res', 105);

[p, n] = findFILE('.', '*_aif_record_new.mat');

load(p{1})

mini_myo_ori = CASignal.m;
v = mean(mini_myo_ori(1:5));
mini_myo = mini_myo_ori - v;

mask_file = 'sr_mask.mat';
mini_im_signal_all = roi_timeseries(squeeze(mini_im(:,:,1,:)), mask_file, 1, 1);
mini_im_signal = mini_im_signal_all.m(4:end) ./ mini_im_signal_all.m(1)
v = mean(mini_im_signal(1:5));
mini_im_signal = mini_im_signal - v;
% ---------------------------------

cd(fullfile(home, full))

full_aif_im = analyze75read('./DebugOutput/aif_upsampled.hdr');
full_aif_signal = analyze75read('./DebugOutput/aif_cin_all_echo0_signal_after_R2StarCorrection_divide_PD_echo0.hdr');

[p, n] = findFILE('.', '*_aif_record_new.mat');

load(p{1})

full_aif = aif_cin_all_echo0_LUTCorrection;

h = figure
hold on
plot(full_aif)
plot(mini_myo*10, 'r');
hold off
title('full aif and mini myo in Gd');
legend('full aif', '10 * mini myo');

cd(home)
saveas(h, 'full_aif_vs_10times_mini_myo_in_Gd', 'fig');

h = figure
hold on
plot(full_aif_signal)
plot(mini_im_signal*10, 'r');
hold off
title('full aif and mini myo in Intensity');
legend('full aif', '10 * mini myo');

cd(home)
saveas(h, 'full_aif_vs_10times_mini_myo_in_IntensitySCC', 'fig');

% -----------------------------------------------------------
% mini, big bolus differences

cd D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150623\stress\meas_MID00039_FID60373_mini__0075_stress_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res

cd D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150601\stress\meas_MID00205_FID53889_MINI__0075_STRESS_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\slep_flow_res

cd D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150506\stress\meas_MID00299_FID45789_MINI_STRESS__0075_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res

cd D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150518\stress\meas_MID00394_FID48382_mini__0075_stress_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res

cd D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150601\stress\meas_MID00205_FID53889_MINI__0075_STRESS_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res

aif_moco_mini = readGTPlusExportImageSeries_Squeeze('.', 101);
perf_moco_mini = readGTPlusExportImageSeries_Squeeze('.', 105);
fig_mini = readGTPlusExportImageSeries_Squeeze('.', 117);
fig_mini = permute(fig_mini, [2 1 3]);

cd D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150623\stress\meas_MID00040_FID60374_big__075_stress_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res

cd D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150601\stress\meas_MID00206_FID53890_BIG__075_STRESS_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\slep_flow_res

cd D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150506\stress\meas_MID00300_FID45790_BIG_STRESS__075_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res

cd D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150518\stress\meas_MID00395_FID48383_big__075_stress_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res

cd D:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150601\stress\meas_MID00206_FID53890_BIG__075_STRESS_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res

aif_moco = readGTPlusExportImageSeries_Squeeze('.', 101);
perf_moco = readGTPlusExportImageSeries_Squeeze('.', 105);
fig = readGTPlusExportImageSeries_Squeeze('.', 117);
fig = permute(fig, [2 1 3]);

PD = 3;

N = size(aif_moco_mini, 4)

N = 40
v = cat(5, aif_moco_mini(:,:,1,PD+1:N), aif_moco(:,:,1,PD+1:N));
v = squeeze(v);
size(v)
figure; imagescn(v, [], [2 N-PD]);

v = cat(5, perf_moco_mini(:,:,1,PD+1:N), perf_moco(:,:,1,PD+1:N));
v = squeeze(v);
size(v)
figure; imagescn(v, [], [2 N-PD]);

figure; imagescn(fig_mini(:,:,2));
figure; imagescn(fig(:,:,2));

% ---------------------------------------------------------------
% SUB, SSFP, linear, nonlinear

cd E:\gtuser\gt_windows_setup\ut\MOCO_EVALUATION\SUB\SUB_SSFP\20150514_09h55m01s_80926

cd E:\gtuser\gt_windows_setup\ut\MOCO_EVALUATION\SUB\SUB_WIP\20150519_08h43m17s_81639

grappa_moco_norm = readGTPlusExportImageSeries_Squeeze('./grappa_flow_res', 105);
grappa_map = readGTPlusExportImageSeries_Squeeze('./grappa_flow_res', 115);

nl_moco_norm = readGTPlusExportImageSeries_Squeeze('./slep_flow_res', 111);
nl_map = readGTPlusExportImageSeries_Squeeze('./slep_flow_res', 115);

figure; imagescn(grappa_moco_norm(:,:,1,:), [], [], [], 4);

figure; imagescn(cat(4, grappa_map, nl_map), [], [2 3]); PerfColorMap;

cd D:\gtuser\gt_windows_setup\ut\MOCO_EVALUATION\SUB\SUB_SSFP\20150514_09h55m01s_80926


cd E:\gtuser\mrprogs\install\DebugOutput_FlowMapping

aif = analyze75read('aif_moco.hdr');
figure; imagescn(aif, [], [], [], 3);
 
perf = analyze75read('perf_moco_0.hdr');
figure; imagescn(perf, [], [], [], 3);

cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150818_A03968\stress\meas_MID00048_FID69679_BIG__075_STRESS_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res

aif = readGTPlusExportImageSeries_Squeeze('.', 1104, 0, 1);

aif1 = aif(:,:,1,1:2:120);
aif1 = squeeze(aif1);
figure; imagescn(aif1, [], [], [], 3);

%% ---------------------------------------
% KLT moco

cd E:\gtuser\mrprogs\install\DebugOutput

perf = analyze75read('Perf_row2.hdr');

KLT0 = analyze75read('KLT_model_Perf_Iter0_row2.hdr');
moco0 = analyze75read('KLT_Perf_moco_after_check_Iter0_row2.hdr');

KLT1 = analyze75read('KLT_model_Perf_Iter1_row2.hdr');
moco1 = analyze75read('KLT_Perf_moco_after_check_Iter1_row2.hdr');

KLT2 = analyze75read('KLT_model_Perf_Iter2_row2.hdr');
moco2 = analyze75read('KLT_Perf_moco_after_check_Iter2_row2.hdr');

perf = flipdim(flipdim(perf, 1), 2);
figure; imagescn(perf, [0 328], [], [], 3);
[I, xlim, ylim, clim, position]=getimagedata(gcf);

KLT0 = flipdim(flipdim(KLT0, 1), 2);
figure; imagescn(KLT0, [0 328], [], [], 3);

KLT1 = flipdim(flipdim(KLT1, 1), 2);
figure; imagescn(KLT1, [0 328], [], [], 3);

KLT2 = flipdim(flipdim(KLT2, 1), 2);
figure; imagescn(KLT2, [0 328], [], [], 3);

moco0 = flipdim(flipdim(moco0, 1), 2);
figure; imagescn(moco0, [0 328], [], [], 3);

moco1 = flipdim(flipdim(moco1, 1), 2);
figure; imagescn(moco1, [0 328], [], [], 3);

moco2 = flipdim(flipdim(moco2, 1), 2);
figure; imagescn(moco2, [0 328], [], [], 3);

setimagescn( xlim, ylim, clim);

save KLT_Moco_Example perf KLT0 KLT1 KLT2 moco0 moco1 moco2

%% ----------------------------

% load flow map

cd E:\gtuser\gt_windows_setup\ut\DualBolus\SSFP\20150901\stress\meas_MID00239_FID76515_BIG__05_STRESS_SSFP_PERF_TPAT3_PF3_4_192x111\grappa_flow_res
a = readGTPlusExportImageSeries_Squeeze('.', 115);

cd E:\gtuser\gt_windows_setup\ut\DualBolus\SSFP\20150901\rest\meas_MID00271_FID76547_BIG__05_REST_SSFP_PERF_TPAT3_PF3_4_192x111\grappa_flow_res
b = readGTPlusExportImageSeries_Squeeze('.', 115);

figure; imagescn(cat(3, a, b)/100, [0 6]); PerfColorMap;

%% -----------------------------

% load AIF figures

cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150901_A04007\meas_MID00456_FID76732_BIG__075_REST_2CC_SEC_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res_new2
a = readGTPlusExportImageSeries_Squeeze('.', 120);
fa = readGTPlusExportImageSeries_Squeeze('.', 115);

cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150901_A04007\meas_MID00457_FID76733_BIG__075_REST_4CC_SEC_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res_new2
b = readGTPlusExportImageSeries_Squeeze('.', 120);
fb = readGTPlusExportImageSeries_Squeeze('.', 115);
c = readGTPlusExportImageSeries_Squeeze('.', 105);


% 2CC
cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150901_A04007\meas_MID00456_FID76732_BIG__075_REST_2CC_SEC_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res2
cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150903_A04012\meas_MID00120_FID77457_BIG__075_2CC_SEC_PERF_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res
cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150923_A04047\meas_MID00229_FID81093_REST_BIG_@2CC_SEC_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res_new2
a = readGTPlusExportImageSeries_Squeeze('.', 120);
fa = readGTPlusExportImageSeries_Squeeze('.', 115);

% 4CC
cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150901_A04007\meas_MID00457_FID76733_BIG__075_REST_4CC_SEC_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res2
cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150923_A04047\meas_MID00230_FID81094_REST_BIG_@4CC_SEC_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res_new2
b = readGTPlusExportImageSeries_Squeeze('.', 120);
fb = readGTPlusExportImageSeries_Squeeze('.', 115);

a = permute(a, [2 1 3]);
b = permute(b, [2 1 3]);

figure; imagescn(a, [], [], 25)
figure; imagescn(b, [], [], 25)


figure; imagescn(fa/100, [0 6], [4 3]); PerfColorMap;

figure; imagescn(cat(3, fa(:,:,:,3), fb(:,:,:,3))/100, [0 6]); PerfColorMap;

%% ------------------------------

% test KLT based aif signal

cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150901_A04007\meas_MID00456_FID76732_BIG__075_REST_2CC_SEC_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150908_A04017\rest\meas_MID00195_FID77924_BIG__075_REST_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res

mask = analyze75read('aif_LV_mask.hdr');
aif_mask_final = analyze75read('AifLVMask_after_Picking.hdr');
aif = analyze75read('input_aif_moco_first_echo_0.hdr');
aifMaskIntensity = analyze75read('aifMaskIntensity.hdr');

ind_all = find(mask(:)>0);
ind = find(aif_mask_final(:)>0);

s = size(aif);

signal_all = zeros(numel(ind_all), s(3));
signal = zeros(numel(ind), s(3));
for p=1:s(3)
    t = aif(:,:,p);
    vt = t(ind(:));
    signal(:, p) = vt;
    
    vt = t(ind_all(:));
    signal_all(:, p) = vt;
end

s = size(signal);
b = signal';
b_std = std(b,0,1);
A = (b'*b)/(s(1)); %clear b;
%clock,
[V,D] = (eig(A)); % This is to diagonalize the covariance matrix
V = (V);
D = (D);
diag(D)

M1 = b*V(:, end);

s = size(aif);
t = 0.5*[0:s(3)-1]
figure
hold on
plot(t, signal_all');
plot(t, b, 'k');
hh = plot(t, M1/sqrt(numel(ind)), 'r--', 'LineWidth', 4);
hold off
box on
legend(hh, 'AIF signal');
title('AIF signal in mask')


cd E:\gtuser\mrprogs\install\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\SUB_Perf\20150831_rest\meas_MID00071_FID100807_PK_rest_Perf_SSFP_GT_TPAT3_PFon_FS_off\grappa_flow_res\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150901_A04007\meas_MID00457_FID76733_BIG__075_REST_4CC_SEC_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150903_A04012\meas_MID00120_FID77457_BIG__075_2CC_SEC_PERF_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150908_A04017\rest\meas_MID00195_FID77924_BIG__075_REST_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\DualBolus\SSFP\20150901_A04005_4mmol_l\stress\meas_MID00239_FID76515_BIG__05_STRESS_SSFP_PERF_TPAT3_PF3_4_192x111\grappa_flow_res\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150908_A04017\stress\meas_MID00175_FID77904_BIG__075_STRESS_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\SUB_Perf\20150831_rest\meas_MID00071_FID100807_PK_rest_Perf_SSFP_GT_TPAT3_PFon_FS_off\grappa_flow_res\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150916_A04034\meas_MID00181_FID79218_BIG__075@2CC_SEC_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res_new\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150601_A03855\rest\meas_MID00227_FID53911_BIG__075_REST_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res_new\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150601_A03855\stress\meas_MID00206_FID53890_BIG__075_STRESS_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\slep_flow_res\DebugOutput_Flow
cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150601_A03855\stress\meas_MID00206_FID53890_BIG__075_STRESS_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res_new\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150601_A03855\rest\meas_MID00227_FID53911_BIG__075_REST_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res_new\DebugOutput

cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150611_A03873\stress\meas_MID00345_FID57138_M_I_N_I__0075_STRESS_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res_new\DebugOutput
cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150629_A03908\stress\meas_MID00128_FID61789_big__075_stress_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res_new\DebugOutput

% low AIF PD, high conc
cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150629_A03908\stress\meas_MID00128_FID61789_big__075_stress_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res_new\DebugOutput

% two injection rate
cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150917_A04037\meas_MID00090_FID79589_BIG__075@4CC_SEC_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res_new\DebugOutput

cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150903_A04012\meas_MID00120_FID77457_BIG__075_2CC_SEC_PERF_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res_new\DebugOutput

cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150916_A04034\meas_MID00180_FID79217_BIG__075@4CC_SEC_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res_new\DebugOutput

cd E:\gtuser\gt_windows_setup\ut\DualBolus\TwoInjectionRate\20150903_A04012\meas_MID00120_FID77457_BIG__075_2CC_SEC_PERF_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res_new2\DebugOutput

aif_0 = analyze75read('aif_for_TwoEcho_T2StartCorrection_0.hdr');
aif_cin_Gd = analyze75read('aif_cin_all_echo0_LUTCorrection.hdr');
aif_cin_all_echo0_signal = analyze75read('aif_cin_all_echo0_signal.hdr');
aif_cin_all_echo1_signal = analyze75read('aif_cin_all_echo1_signal.hdr');
aif_cin_all_echo0_signal_after_R2StarCorrection = analyze75read('aif_cin_all_echo0_signal_after_R2StarCorrection.hdr');
aif_cin_all_echo0_OverPD_after_R2StarCorrection = analyze75read('aif_cin_all_echo0_OverPD_after_R2StarCorrection.hdr');
aif_LUT = analyze75read('aif_cin_LUT_flash_pd_flash_sr.hdr');
aif_cin_all_R2Star = analyze75read('aif_cin_all_R2Star.hdr');
aif_cin_all_R2Star_SLEP = analyze75read('aif_cin_all_R2Star_SLEP.hdr');
aif_PD = analyze75read('aifPD_for_TwoEcho_T2StartCorrection_0.hdr');
aif_mask = analyze75read('aif_LV_mask_for_TwoEcho_T2StartCorrection_0.hdr');
aif_mask_final = analyze75read('AifLVMask_after_Picking.hdr');
PDMaskIntensity = analyze75read('PDMaskIntensity.hdr');
pdPicked = analyze75read('pdPicked.hdr');
aifPicked = analyze75read('aifPicked.hdr');
aifMaskIntensity = analyze75read('aifMaskIntensity.hdr');
aif0 = analyze75read('AIF_input_for_moco_0_MAG.hdr');
aif_LUT = analyze75read('aif_cin_LUT_flash_pd_flash_sr.hdr');


aif_PD_echo0 = analyze75read('input_aif_PD_0.hdr');
aif_PD_echo1 = analyze75read('input_aif_PD_1.hdr');

figure; imagescn(cat(4, aif_PD_echo0, aif_PD_echo1), [], [2 3]);
figure; imagescn(aif_0, [], [], [], 3);
figure; imagescn(aif_mask_final, [], [], [], 3);

ind = find(aif_mask_final(:)>0);

PD = aif_PD_echo0(:, :, 1);
PD_v = PD(ind(:));
mean(PD_v)

aif_s0 = zeros(size(aif_0, 3), 1);

for ii=1:size(aif_0, 3)
    v = aif_0(:,:,ii);
    p = v(ind(:));
    aif_s0(ii) = mean(p);
end

figure;
hold on
plot(aif_cin_all_echo0_signal);
plot(aif_cin_all_echo1_signal, 'r');
plot(aif_cin_all_echo0_signal_after_R2StarCorrection, 'k');
plot(aif_s0, 'g-');
hold off

figure;
hold on
plot(aif_cin_all_R2Star);
plot(aif_cin_all_R2Star_SLEP, 'r');
hold off

figure;
hold on
plot(aif_cin_all_echo0_OverPD_after_R2StarCorrection)
plot(aif_cin_all_echo0_signal_after_R2StarCorrection/mean(PD_v), 'r-.')
hold off

figure;
hold on
plot(aif_cin_Gd)
hold off

figure
plot([0:0.1:20], aif_LUT);
xlabel('Gd')
ylabel('SR/PD')

figure;
hold on
plot(aifMaskIntensity', 'k', 'LineWidth', 1);
plot(abs(aif_cin_all_echo0_signal), 'r--', 'LineWidth', 4);
plot(aifPicked', 'r');
plot(mean(aifPicked', 2), 'g--', 'LineWidth', 4);
hold off





r2s = log(0.6/0.4726) / (1.65-1.16);

0.6*exp(1.16*r2s)

a1 = 2863
a2 = 2853
te0 = 0.6
te1 = 1.49
r2s = log(a1/a2) / (te1-te0)
a1*exp(te0*r2s)




aif0_s = aif0(:,:, 4:end);
aif0_s2 = reshape(aif0_s, [64*48 57]);

s = size(aif0_s2);
b = aif0_s2;
b_std = std(b,0,1);
A = (b'*b)/(s(1)); %clear b;
[V,D] = (eig(A)); % This is to diagonalize the covariance matrix
V = (V);
D = (D);
diag(D)

M1 = b*V(:, end);

figure
hold on
plot(b);
plot(M1/sqrt(42), 'r--', 'LineWidth', 4);
hold off



cd E:\gtuser\gt_windows_setup\ut\SUB_Perf\20150831_rest\meas_MID00071_FID100807_PK_rest_Perf_SSFP_GT_TPAT3_PFon_FS_off\grappa_flow_res
cd E:\gtuser\gt_windows_setup\ut\SUB_Perf\20150423\20150423_10h09m24s_77558\grappa_flow_res

a = readGTPlusExportImageSeries_Squeeze('.', 115);
figure; imagescn(a/100, [0 6], []); PerfColorMap;


a0 = analyze75read('CASignal_Perf_0.hdr');
a1 = analyze75read('CASignal_Perf_1.hdr');
a2 = analyze75read('CASignal_Perf_2.hdr');

figure; imagescn(cat(4, a0, a1, a2), [], [], [], 3);


%% new vs. old

cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150908_A04017\stress\meas_MID00175_FID77904_BIG__075_STRESS_FLASH_PERF_TPAT3_PF3_4_192x111\grappa_flow_res\DebugOutput

cd E:\gtuser\gt_windows_setup\ut\DualBolus\SSFP\20150911_A04027_4mmol_l\stress\meas_MID00047_FID78472_BIG__05_STRESS_SSFP_PERF_TPAT3_PF3_4_192x111\grappa_flow_res\DebugOutput

cd E:\gtuser\gt_windows_setup\ut\DualBolus\SSFP\20150901_A04005_4mmol_l\stress\meas_MID00239_FID76515_BIG__05_STRESS_SSFP_PERF_TPAT3_PF3_4_192x111\grappa_flow_res\DebugOutput

% 1.5T
cd E:\gtuser\gt_windows_setup\ut\SUB_Perf\20150831_rest\meas_MID00071_FID100807_PK_rest_Perf_SSFP_GT_TPAT3_PFon_FS_off\grappa_flow_res\DebugOutput

cd E:\gtuser\gt_windows_setup\ut\SUB_Perf\cases\E07497\rest\20150818_09h10m39s_97598\grappa_flow_res\DebugOutput

SR_Norm_new = analyze75read('SRNorm_0.hdr');
LUT_new = analyze75read('Perf_T1_Correction_LUT_perf_ssfp_PD_flash.hdr');
aif_LUT_new = analyze75read('aif_cin_LUT_flash_pd_flash_sr.hdr');
PD_new = analyze75read('input_perf_PD_0.hdr');

cd E:\gtuser\gt_windows_setup\ut\DualBolus\FLASH\20150728_A03964\stress\meas_MID00124_FID68823_BIG__075_STRESS_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\grappa_flow_res\DebugOutput

cd E:\gtuser\gt_windows_setup\ut\DualBolus\SSFP\20150731_A03965\stress\meas_MID00038_FID69056_BIG__05_STRESS_SSFP_PERF_TPAT3_PF3_4_192x111\grappa_flow_res\DebugOutput

% 1.5T
cd E:\gtuser\gt_windows_setup\ut\SUB_Perf\20150423\20150423_10h09m24s_77558\grappa_flow_res\DebugOutput

SR_Norm_old = analyze75read('SRNorm_0.hdr');
LUT_old = analyze75read('Perf_T1_Correction_LUT_perf_ssfp_PD_flash.hdr');
aif_LUT_old = analyze75read('aif_cin_LUT_flash_pd_flash_sr.hdr');
PD_old = analyze75read('input_perf_PD_0.hdr');

figure; imagescn(cat(4, SR_Norm_old, SR_Norm_new), [], [], [], 3);

figure; imagescn(cat(4, PD_old, PD_new), [], [], [], 3);

save SR_Norm_Old_New SR_Norm_new SR_Norm_old

figure
hold on
plot(0:0.1:20, aif_LUT_new)
plot(0:0.1:20, aif_LUT_old, 'r')
hold off
legend('new', 'old')

figure
hold on
plot(0:0.01:20, LUT_new)
plot(0:0.01:20, LUT_old, 'r')
hold off
legend('new', 'old')
xlabel('Gd')
ylabel('SR/PD')

%% viewing

for ii=1:num  
   
    if ( ~isempty(strfind(subdir{ii}, 'mini')) | ~isempty(strfind(subdir{ii}, 'MINI'))  )
        continue;
    end
    
    subdir{ii}  

    cd(fullfile(home, subdir{ii}))
  
    try
        [flowmaps, acq_time, physio_time] = readGTPlusExportImageSeries('./grappa_flow_res_new3', 115, 1);
        flowmaps = squeeze(flowmaps/100);
        h = figure; imagescn(flowmaps, [0 6], [4 3]); PerfColorMap;
        saveas(h, 'flowmaps', 'fig');
    catch
    end
    
    try
        [figs, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze('./grappa_flow_res_new3', 120, 1);
        figs = permute(squeeze(figs), [2 1 3]);
        h = figure; imagescn(figs, [], [], 20);
        saveas(h, 'GdFigs', 'fig');
    catch
    end
end

%% remove folder

for pp=1:numel(home_all)
    
    home = home_all{pp}
    
%     [subdir, num] = FindAllEndDirectoryExclusive(home, exDir)
    
    home_ut = home(length(UTDir)+2:end)
    UTCase_Flow{1} = home_ut;

    % ------------------------------------------------------------

    exDir = {'ICERecon', 'grappa', 'TXMapping', 'moco', 'mocoSyn', 'mocoPS' , 'mocoPSSyn', 'seg', 'magFitting', 'slep_res', 'slep_cloud_res', 'grappa_res', 'ICE', 'DebugOutput', 'PhantomT1maps', 'slep_flow_res', 'grappa_flow_res', 'slep_cloud_flow_res', 'slep_cloud_flow_res2', 'T2Map'};
    [subdir, num] = FindAllEndDirectoryExclusive(home, exDir)

    % ------------------------------------------------------------
    
    startI = 1;
%     if(pp==3)
%         startI = 18;
%     end
    
    for ii=startI:num  

        if ( ~isempty(strfind(subdir{ii}, 'mini')) | ~isempty(strfind(subdir{ii}, 'MINI')) | ~isempty(strfind(subdir{ii}, 'Mini')) )
%             continue;
        end
        
        fullPath = fullfile(home, subdir{ii}, 'grappa_flow_res_BTEX20_2')
        try
            rmdir(fullPath, 's');
        catch
        end
        
        fullPath = fullfile(home, subdir{ii}, 'grappa_flow_res_BTEX20_3')
        try
            rmdir(fullPath, 's');
        catch
        end
        
        fullPath = fullfile(home, subdir{ii}, 'grappa_flow_res_BTEX20_4')
        try
            rmdir(fullPath, 's');
        catch
        end
        
        fullPath = fullfile(home, subdir{ii}, 'grappa_flow_res_new')
        try
            rmdir(fullPath, 's');
        catch
        end
        
        fullPath = fullfile(home, subdir{ii}, 'grappa_flow_res_new_2')
        try
            rmdir(fullPath, 's');
        catch
        end
        
        fullPath = fullfile(home, subdir{ii}, 'grappa_flow_res_new2')
        try
            rmdir(fullPath, 's');
        catch
        end
        
        fullPath = fullfile(home, subdir{ii}, 'grappa_flow_res_BTEX20_8')
        try
            rmdir(fullPath, 's');
        catch
        end 
        
        fullPath = fullfile(home, subdir{ii}, 'grappa_flow_res_new3')
        try
            rmdir(fullPath, 's');
        catch
        end
        
        fullPath = fullfile(home, subdir{ii}, 'grappa_flow_res_PSIR')
        try
            rmdir(fullPath, 's');
        catch
        end
        
        fullPath = fullfile(home, subdir{ii}, 'grappa_flow_res_PSIR_correctedAIF')
        try
            rmdir(fullPath, 's');
        catch
        end
        
        fullPath = fullfile(home, subdir{ii}, 'grappa_flow_res_BTEX20')
        try
            rmdir(fullPath, 's');
        catch
        end
    end
end
