

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

    for ii=startI:num  

        if ( ~isempty(strfind( lower(subdir{ii}), 'mini')) | ~isempty(strfind(subdir{ii}, 'MINI')) | ~isempty(strfind(subdir{ii}, 'Mini')) )
            continue;
        end
        
        if ( ~isempty(strfind( lower(subdir{ii}), 't1map')) )
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

        [names_dat, num_dat] = findFILE(fullfile(home, subdir{ii}), '*.dat')
        
        if(num_dat==0)
            [names_dat, num_dat] = findFILE(fullfile(home, subdir{ii}), '*.h5')
        end
        
        if(num_dat==0)
            continue;
        end
        
        try
            rmdir(fullfile(home, subdir{ii}, 'Dotarem_r1_3p94_r2_4p73'), 's');
        catch
        end
        
        try
            rmdir(fullfile(home, subdir{ii}, 'grappa_flow_res_BTEX20_TwoCompWithShifts_Dotarem_r1_3p94_r2_4p73'), 's');
        catch
        end
        
        try
            rmdir(fullfile(home, subdir{ii}, 'grappa_flow_res'), 's');
        catch
        end

        try
            rmdir(fullfile(home, subdir{ii}, 'grappa_flow_res_BTEX20_TwoCompWithShifts_r1_3p94_r2_4p73'), 's');
        catch
        end
        
        try
            rmdir(fullfile(home, subdir{ii}, 'grappa_flow_res_BTEX20_TwoCompWithShifts_Dotarem_r1_3p94_r2_4p73'), 's');
        catch
        end
        
        try
            rmdir(fullfile(home, subdir{ii}, 'Dotarem_r1_r2_increased'), 's');
        catch
        end
        
        try
            rmdir(fullfile(home, subdir{ii}, 'Dotarem_r1_r2_reduced'), 's');
        catch
        end
        
        delete(fullfile(home, subdir{ii}, 'slep_flow_res_BTEX20_5', 'DebugOutput/data_*'));
        delete(fullfile(home, subdir{ii}, 'slep_flow_res_BTEX20_5', 'DebugOutput/workOrder_ref*'));
        delete(fullfile(home, subdir{ii}, 'slep_flow_res_BTEX20_5', 'DebugOutput/ref_*'));
        delete(fullfile(home, subdir{ii}, 'slep_flow_res_BTEX20_5', 'DebugOutput/incomingKSpace_*'));
        delete(fullfile(home, subdir{ii}, 'slep_flow_res_BTEX20_5', 'DebugOutput/incomingRef_*'));
        
        delete(fullfile(home, subdir{ii}, 'slep_flow_res_BTEX20_TwoCompWithShifts', 'DebugOutput/data_*'));
        delete(fullfile(home, subdir{ii}, 'slep_flow_res_BTEX20_TwoCompWithShifts', 'DebugOutput/workOrder_ref*'));
        delete(fullfile(home, subdir{ii}, 'slep_flow_res_BTEX20_TwoCompWithShifts', 'DebugOutput/workOrder_data*'));
        delete(fullfile(home, subdir{ii}, 'slep_flow_res_BTEX20_TwoCompWithShifts', 'DebugOutput/ref_*'));
        delete(fullfile(home, subdir{ii}, 'slep_flow_res_BTEX20_TwoCompWithShifts', 'DebugOutput/incomingKSpace_*'));
        delete(fullfile(home, subdir{ii}, 'slep_flow_res_BTEX20_TwoCompWithShifts', 'DebugOutput/incomingRef_*'));
    end

end

