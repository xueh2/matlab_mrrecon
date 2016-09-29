function ReviewPerfusionRecon_OneFolder(UT_Dir, home, res_dir, seriesNo)
% ReviewPerfusionRecon_OneFolder('F', 'DualBolus\FFR', 'test_moco', 120)

%% 
set_UT_Dir(UT_Dir)
UTDir = getenv('GTPLUS_UT_DIR')

home_all = fullfile(UTDir, home)

exDir = {'ICERecon', 'grappa', 'TXMapping', 'moco', 'mocoSyn', 'mocoPS' , 'mocoPSSyn', 'seg', 'magFitting', 'slep_res', 'slep_cloud_res', 'grappa_res', 'ICE', 'DebugOutput', 'PhantomT1maps', 'slep_flow_res', 'grappa_flow_res', 'slep_cloud_flow_res', 'slep_cloud_flow_res2', 'T2Map', 'test', 'Dotarem_r1_r2_reduced', 'grappa_flow_res_BTEX20_TwoCompWithShifts', 'slep_flow_res_BTEX20_TwoCompWithShifts'};
[subdir, num] = FindAllEndDirectoryExclusive(home_all, exDir)

startI = 1;

for ii=startI:num  

    if ( ~isempty(strfind( lower(subdir{ii}), 'mini')) | ~isempty(strfind(subdir{ii}, 'MINI')) | ~isempty(strfind(subdir{ii}, 'Mini')) )
        continue;
    end

    if ( ~isempty(strfind( lower(subdir{ii}), 't1map')) )
        continue;
    end

    resDir = fullfile(home_all, subdir{ii}, res_dir)
    
    try
        a = readGTPlusExportImageSeries_Squeeze(resDir, seriesNo);
    catch
        closeall;
        continue;
    end
    
    if(seriesNo==120)
        a = permute(a, [2 1 3]);
        figure; imagescn(a, [], [], 24);
    else    
        n = numel(size(a));
        figure; imagescn(a, [], [], [], n);
    end
    
    pause;
    closeall
end
