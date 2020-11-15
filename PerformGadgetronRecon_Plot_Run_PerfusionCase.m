
function PerformGadgetronRecon_Plot_Run_PerfusionCase(perf_cases, rest_cases, dataDir, resDir, host, configNamePreset, checkProcessed, flow_windowing, only_reviewing, only_processing, only_making_figures, startRemoteGT, delete_old_res, gt_port)
% PerformGadgetronRecon_Plot_Run_PerfusionCase(perf_cases, rest_cases, dataDir, resDir, host, configNamePreset, checkProcessed, flow_windowing, only_reviewing, only_processing, only_making_figures, startRemoteGT, delete_old_res, gt_port)

if(nargin < 8)
    flow_windowing = [0 8];
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

if(nargin < 14)
    gt_port = [];
end

if(strcmp(host, 'beast'))
    host = '137.187.135.157';
end

if(strcmp(host, 'gt1'))
    host = '137.187.135.169';
end

if(strcmp(host, 'gt2'))
    host = '137.187.135.238';
end

noise_dat_processed = [];

startN = 1
endN = size(perf_cases, 1)
for ii=startN:endN
    
    perf_xml = find_perf_xml(configNamePreset, perf_cases{ii, 2});
    disp([num2str(ii-startN+1) ' out of ' num2str(endN-startN+1) ' - ' perf_cases(ii, 2:3) ' - ' perf_xml]);
    
    if(~only_making_figures)
        if(~only_processing)
            try
                [h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir, perf_cases{ii, 2}, perf_cases{ii, 3}, flow_windowing, only_reviewing, checkProcessed);
            catch
                h_flow_stress = 0;
                has_stress = 0;
                has_rest = 0;
            end
            
            if(h_flow_stress==0 | ~ishandle(h_flow_stress) | h_flow_rest==0 | ~ishandle(h_flow_rest))
                disp(['--> Processing ... ']);
                
                if(~has_stress)
                    PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, perf_cases(ii, 2), host, resDir, 0, delete_old_res, startRemoteGT, {perf_xml}, noise_dat_processed, gt_port);
                end

                if(~has_rest)
                    PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, perf_cases(ii, 3), host, resDir, 0, delete_old_res, startRemoteGT, {perf_xml}, noise_dat_processed, gt_port);
                end
%                 PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, perf_cases(ii, 2:3), host, resDir, checkProcessed, delete_old_res, startRemoteGT, {perf_xml});

                [h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir, perf_cases{ii, 2}, perf_cases{ii, 3}, flow_windowing, only_reviewing, checkProcessed);
            end
        else
            has_stress = 0;
            has_rest = 0;
            noise_dat_processed = [];
            
            if(~has_stress)
                PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, perf_cases(ii, 2), host, resDir, 0, delete_old_res, startRemoteGT, {perf_xml}, noise_dat_processed, gt_port);
            end

            if(~has_rest)
                PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, perf_cases(ii, 3), host, resDir, 0, delete_old_res, startRemoteGT, {perf_xml}, noise_dat_processed, gt_port);
            end

            %[h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir, perf_cases{ii, 2}, perf_cases{ii, 3}, flow_windowing, only_reviewing, 0);
        end
    end
    
    if(~only_processing)
        try
            [h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir, perf_cases{ii, 2}, perf_cases{ii, 3}, flow_windowing, only_reviewing, checkProcessed);
        catch
        end
        
        if(has_stress<0)
            PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, perf_cases(ii, 2), host, resDir, 0, delete_old_res, startRemoteGT, {perf_xml}, noise_dat_processed, gt_port);
        end

        if(has_rest<0)
            PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, perf_cases(ii, 3), host, resDir, 0, delete_old_res, startRemoteGT, {perf_xml}, noise_dat_processed, gt_port);
        end
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
    
    if(isempty(strfind(caseName, 'Perfusion')))
        caseName = rest_cases{ii, 1};
    end
    
    perf_xml = find_perf_xml(configNamePreset, caseName);
    disp([num2str(ii-startN+1) ' out of ' num2str(endN-startN+1) ' - ' caseName ' - ' perf_xml]);
    
    if(~only_making_figures)
        if(~only_processing)
            try
                h_flow_rest = PerformGadgetronRecon_Plot_PerfusionCase(resDir, caseName, only_reviewing);
                if(h_flow_rest==0)
                    PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, {caseName}, host, resDir, 0, delete_old_res, startRemoteGT, {perf_xml}, noise_dat_processed, gt_port);

                    PerformGadgetronRecon_Plot_PerfusionCase(resDir, caseName, only_reviewing);
                end
            catch
                disp(['--> Processing ... ']);
                PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, {caseName}, host, resDir, 0, delete_old_res, startRemoteGT, {perf_xml}, noise_dat_processed, gt_port);

                PerformGadgetronRecon_Plot_PerfusionCase(resDir, caseName, only_reviewing);
            end
        else
            PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, {caseName}, host, resDir, 0, delete_old_res, startRemoteGT, {perf_xml}, noise_dat_processed, gt_port);

%             PerformGadgetronRecon_Plot_PerfusionCase(resDir, caseName, only_reviewing);
        end
    else
        PerformGadgetronRecon_Plot_PerfusionCase(resDir, caseName, only_reviewing, checkProcessed);
    end
    
    if(only_reviewing)
        pause;
    end
    closeall
    
    disp(['====================================================================================================']);
end

end

function perf_xml = find_perf_xml(configNames, dataName)

    if(~iscell(configNames))
        perf_xml = configNames(1);
    end
    
    R3_xml = zeros(numel(configNames), 1);
    for n=1:numel(configNames)
        ind = strfind(configNames{n}, 'AIFR3');
        if(isempty(ind))
            R3_xml(n) = 0;
        else
            R3_xml(n) = 1;
        end
    end
    
    ind_r3 = find(R3_xml==1);
    ind_r2 = find(R3_xml==0);
    
    ind = strfind(dataName, 'R3');
    if(isempty(ind))
        perf_xml = configNames{ind_r2(1)};
    else
        perf_xml = configNames{ind_r3(1)};
    end
end
