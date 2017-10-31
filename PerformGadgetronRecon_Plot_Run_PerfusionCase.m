
function PerformGadgetronRecon_Plot_Run_PerfusionCase(perf_cases, rest_cases, dataDir, resDir, host, configNamePreset, checkProcessed, flow_windowing, only_reviewing, only_processing, only_making_figures, startRemoteGT, delete_old_res)
% PerformGadgetronRecon_Plot_Run_PerfusionCase(perf_cases, rest_cases, dataDir, resDir, host, configNamePreset, checkProcessed, flow_windowing, only_reviewing, only_processing, only_making_figures, startRemoteGT, delete_old_res)

if(nargin < 8)
    flow_windowing = [0 6];
end

if(nargin < 9)
    only_reviewing = 0;
end

if(nargin < 10)
    only_processing = 0;
end

if(nargin < 11)
    only_making_figures = 0;
end

if(nargin < 12)
    startRemoteGT = 1;
end

if(nargin < 13)
    delete_old_res = 1;
end

startN = 1
endN = size(perf_cases, 1)
for ii=startN:endN
    disp([num2str(ii-startN+1) ' out of ' num2str(endN-startN+1) ' - ' perf_cases(ii, 2:3)]);
    
    if(~only_making_figures)
        if(~only_processing)
            try
                [h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir, perf_cases{ii, 2}, perf_cases{ii, 3}, flow_windowing, only_reviewing);
            catch
                disp(['--> Processing ... ']);
                PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, perf_cases(ii, 2:3), host, resDir, checkProcessed, delete_old_res, startRemoteGT, {configNamePreset});

                [h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir, perf_cases{ii, 2}, perf_cases{ii, 3}, flow_windowing, only_reviewing);
            end
        else
            has_stress = 0;
            has_rest = 0;
        end

        if(~has_stress)
            PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, perf_cases(ii, 2), host, resDir, checkProcessed, delete_old_res, startRemoteGT, {configNamePreset});
        end

        if(~has_rest)
            PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, perf_cases(ii, 3), host, resDir, checkProcessed, delete_old_res, startRemoteGT, {configNamePreset});
        end
    end
    
    try
        [h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir, perf_cases{ii, 2}, perf_cases{ii, 3}, flow_windowing, only_reviewing, checkProcessed);
    catch
    end
    
    if(only_reviewing)
        pause;
    end
    closeall
    
    disp(['====================================================================================================']);
end

startN = 1
endN = size(rest_cases, 1)
for ii=startN:endN
    caseName = rest_cases{ii, 2};
    disp([num2str(ii-startN+1) ' out of ' num2str(endN-startN+1) ' - ' caseName]);
    
    if(~only_making_figures)
        if(~only_processing)
            try
                h_flow_rest = PerformGadgetronRecon_Plot_PerfusionCase(resDir, caseName, only_reviewing);
                if(h_flow_rest==0)
                    PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, {caseName}, host, resDir, checkProcessed, delete_old_res, startRemoteGT, {configNamePreset});

                    PerformGadgetronRecon_Plot_PerfusionCase(resDir, caseName, only_reviewing);
                end
            catch
                disp(['--> Processing ... ']);
                PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, {caseName}, host, resDir, checkProcessed, delete_old_res, startRemoteGT, {configNamePreset});

                PerformGadgetronRecon_Plot_PerfusionCase(resDir, caseName, only_reviewing);
            end
        else
            PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, {caseName}, host, resDir, checkProcessed, delete_old_res, startRemoteGT, {configNamePreset});

            PerformGadgetronRecon_Plot_PerfusionCase(resDir, caseName, only_reviewing);
        end
    else
        PerformGadgetronRecon_Plot_PerfusionCase(resDir, caseName, only_reviewing);
    end
    
    if(only_reviewing)
        pause;
    end
    closeall
    
    disp(['====================================================================================================']);
end
