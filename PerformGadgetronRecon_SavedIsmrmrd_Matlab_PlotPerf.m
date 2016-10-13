
function PerformGadgetronRecon_SavedIsmrmrd_Matlab_PlotPerf(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, fixed_HCT, prefix)
% PerformGadgetronRecon_SavedIsmrmrd_Matlab_PlotPerf(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, fixed_HCT, prefix)

flow_window = [0 6];
scalingFactor = 10;

nV = numel(PerfTable(1, :));

num_column = size(PerfTable, 2);
num = size(PerfTable, 1)-1;
for n=1:num
       
    stressCase = PerfTable{n+1, stress_column};
    restCase = PerfTable{n+1, rest_column};
    if(ischemia_column<=num_column)
        Category = PerfTable{n+1, ischemia_column};
    else
        Category = 1;
    end
    
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]); 
    PerfTable{n+1, :}
    disp(['==================================================================']);  
       
    [configName, scannerID, pID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(stressCase);

    figDir = fullfile(resDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_Figure']) 
    
    roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_ROI'])
        
    if(~isempty(fixed_HCT))
        if(fixed_HCT>=0 & fixed_HCT<=1)
            HCT = fixed_HCT;
        end
    else
        HCT = PerfTable{n+1, hct_column};
        if(isnan(HCT))
            HCT = 0;
        else
            if(HCT>1)
                HCT = HCT/100;
            end
        end
    end
    
    disp([num2str(n) ' out of ' num2str(num) ' - HCT : ' num2str(HCT)]); 
   
    hct_str = num2str(HCT);
    ind = find(hct_str=='.');
    if(~isempty(ind))
        hct_str(ind(:)) = 'p';
    end
    suffix = [study_dates '_hct' hct_str];        
            
    stressMat1 = load(fullfile(resDir, study_dates, stressCase, [prefix '_' suffix '_0.mat']));
    stressMat2 = load(fullfile(resDir, study_dates, stressCase, [prefix '_' suffix '_1.mat']));
    stressMat3 = load(fullfile(resDir, study_dates, stressCase, [prefix '_' suffix '_2.mat']));

    restMat1 = load(fullfile(resDir, study_dates, restCase, [prefix '_' suffix '_0.mat']));
    restMat2 = load(fullfile(resDir, study_dates, restCase, [prefix '_' suffix '_1.mat']));
    restMat3 = load(fullfile(resDir, study_dates, restCase, [prefix '_' suffix '_2.mat']));
                  
    res_stress = PerformGadgetronRecon_Matlab_ROIValues_OneCase(stressMat1, stressMat2, stressMat3, HCT);
    res_rest = PerformGadgetronRecon_Matlab_ROIValues_OneCase(restMat1, restMat2, restMat3, HCT);
    
    if(isFileExist(roiDir))
        cd(roiDir)
    end
    
    slc = size(res_rest.flow, 3);
    
    if(size(res_stress.flow, 2)==size(res_rest.flow,2))
        h = figure('Name','Flow maps','NumberTitle','off'); imagescn(cat(3, res_stress.flow, res_rest.flow), flow_window, [2 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE Visf','NumberTitle','off'); imagescn(cat(3, res_stress.Visf, res_rest.Visf), [0 100], [2 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE PS','NumberTitle','off'); imagescn(cat(3, res_stress.PS, res_rest.PS), [0 10], [2 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE E','NumberTitle','off'); imagescn(cat(3, res_stress.E, res_rest.E), [0 2], [2 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE Vp','NumberTitle','off'); imagescn(cat(3, res_stress.Vp, res_rest.Vp), [0 20], [2 slc], scalingFactor); PerfColorMap;    
    else
        h = figure('Name','Flow maps, stress','NumberTitle','off'); imagescn(res_stress.flow, flow_window, [1 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE Visf, stress','NumberTitle','off'); imagescn(res_stress.Visf, [0 100], [1 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE PS, stress','NumberTitle','off'); imagescn(res_stress.PS, [0 10], [1 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE E, stress','NumberTitle','off'); imagescn(res_stress.E, [0 2], [1 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE Vp, stress','NumberTitle','off'); imagescn(res_stress.Vp, [0 20], [1 slc], scalingFactor); PerfColorMap;
        
        h = figure('Name','Flow maps, rest','NumberTitle','off'); imagescn(res_rest.flow, flow_window, [1 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE Visf, rest','NumberTitle','off'); imagescn(res_rest.Visf, [0 100], [1 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE PS, rest','NumberTitle','off'); imagescn(res_rest.PS, [0 10], [1 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE E, rest','NumberTitle','off'); imagescn(res_rest.E, [0 2], [1 slc], scalingFactor); PerfColorMap;
        h = figure('Name','PDE Vp, rest','NumberTitle','off'); imagescn(res_rest.Vp, [0 20], [1 slc], scalingFactor); PerfColorMap;
    end
    pause;
    closeall
end

end

function res = PerformGadgetronRecon_Matlab_ROIValues_OneCase(a1, a2, a3, HCT)
% res = PerformGadgetronRecon_Matlab_ROIValues_OneCase(a1, a2, a3)
% given the ROIs s1/s2/s3, get the flow and other values
% res has [flow, Ki, E, Visf, Vp, PS, SD, Ki_MF, Ki_Fermi, Ki_TwoCompExp, Ki_BTEX]

if(isfield(a1, 'flowmaps_grappa_PSIR'))
    
    res.flow = cat(3, a1.flowmaps_grappa_PSIR(:,:,end), a2.flowmaps_grappa_PSIR(:,:,end), a3.flowmaps_grappa_PSIR(:,:,end));

    EMap1 = a1.Ki_whole_grappa_PSIR(:,:,4) ./ (a1.flowmaps_grappa_PSIR(:,:,end)+eps);
    EMap2 = a2.Ki_whole_grappa_PSIR(:,:,4) ./ (a2.flowmaps_grappa_PSIR(:,:,end)+eps);
    EMap3 = a3.Ki_whole_grappa_PSIR(:,:,4) ./ (a3.flowmaps_grappa_PSIR(:,:,end)+eps);
    
    res.E = cat(3, EMap1, EMap2, EMap3);
    res.Visf = cat(3, a1.grappa_interVolumeMap_grappa_PSIR(:,:,end), a2.grappa_interVolumeMap_grappa_PSIR(:,:,end), a3.grappa_interVolumeMap_grappa_PSIR(:,:,end));

    res.Vp = cat(3, a1.blood_volume_maps_grappa_PSIR(:,:,end), a2.blood_volume_maps_grappa_PSIR(:,:,end), a3.blood_volume_maps_grappa_PSIR(:,:,end)) * (1-HCT);
    res.PS = cat(3, a1.PS_maps_grappa_PSIR(:,:,end), a2.PS_maps_grappa_PSIR(:,:,end), a3.PS_maps_grappa_PSIR(:,:,end));
    res.SD = cat(3, a1.SD_maps_grappa_PSIR(:,:,end), a2.SD_maps_grappa_PSIR(:,:,end), a3.SD_maps_grappa_PSIR(:,:,end));
    res.Ki = cat(4, a1.Ki_whole_grappa_PSIR, a2.Ki_whole_grappa_PSIR, a3.Ki_whole_grappa_PSIR);

elseif(isfield(a1, 'flowmaps_grappa_PSIR_OnlyGlobalSearch'))
    
    res.flow = cat(3, a1.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end), a2.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end), a3.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end));

    EMap1 = a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4) ./ (a1.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end)+eps);
    EMap2 = a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4) ./ (a2.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end)+eps);
    EMap3 = a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4) ./ (a3.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end)+eps);
    
    res.E = cat(3, EMap1, EMap2, EMap3);
    res.Visf = cat(3, a1.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch(:,:,end), a2.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch(:,:,end), a3.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch(:,:,end));

    res.Vp = cat(3, a1.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch(:,:,end), a2.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch(:,:,end), a3.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch(:,:,end)) * (1-HCT);
    res.PS = cat(3, a1.PS_maps_grappa_PSIR_OnlyGlobalSearch(:,:,end), a2.PS_maps_grappa_PSIR_OnlyGlobalSearch(:,:,end), a3.PS_maps_grappa_PSIR_OnlyGlobalSearch(:,:,end));
    res.SD = cat(3, a1.SD_maps_grappa_PSIR_OnlyGlobalSearch(:,:,end), a2.SD_maps_grappa_PSIR_OnlyGlobalSearch(:,:,end), a3.SD_maps_grappa_PSIR_OnlyGlobalSearch(:,:,end));
    res.Ki = cat(4, a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch, a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch, a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch);

elseif(isfield(a1, 'flowmaps_grappa_PSIR_without_R2Star'))
    
    res.flow = cat(3, a1.flowmaps_grappa_PSIR_without_R2Star(:,:,end), a2.flowmaps_grappa_PSIR_without_R2Star(:,:,end), a3.flowmaps_grappa_PSIR_without_R2Star(:,:,end));

    EMap1 = a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,4) ./ (a1.flowmaps_grappa_PSIR_without_R2Star(:,:,end)+eps);
    EMap2 = a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,4) ./ (a2.flowmaps_grappa_PSIR_without_R2Star(:,:,end)+eps);
    EMap3 = a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,4) ./ (a3.flowmaps_grappa_PSIR_without_R2Star(:,:,end)+eps);
    
    res.E = cat(3, EMap1, EMap2, EMap3);
    res.Visf = cat(3, a1.grappa_interVolumeMap_grappa_PSIR_without_R2Star(:,:,end), a2.grappa_interVolumeMap_grappa_PSIR_without_R2Star(:,:,end), a3.grappa_interVolumeMap_grappa_PSIR_without_R2Star(:,:,end));

    res.Vp = cat(3, a1.blood_volume_maps_grappa_PSIR_without_R2Star(:,:,end), a2.blood_volume_maps_grappa_PSIR_without_R2Star(:,:,end), a3.blood_volume_maps_grappa_PSIR_without_R2Star(:,:,end)) * (1-HCT);
    res.PS = cat(3, a1.PS_maps_grappa_PSIR_without_R2Star(:,:,end), a2.PS_maps_grappa_PSIR_without_R2Star(:,:,end), a3.PS_maps_grappa_PSIR_without_R2Star(:,:,end));
    res.SD = cat(3, a1.SD_maps_grappa_PSIR_without_R2Star(:,:,end), a2.SD_maps_grappa_PSIR_without_R2Star(:,:,end), a3.SD_maps_grappa_PSIR_without_R2Star(:,:,end));
    res.Ki = cat(4, a1.Ki_whole_grappa_PSIR_without_R2Star, a2.Ki_whole_grappa_PSIR_without_R2Star, a3.Ki_whole_grappa_PSIR_without_R2Star);

end

res.flow = flipdim(flipdim(res.flow,1), 2);
res.E = flipdim(flipdim(res.E,1), 2);
res.Visf = flipdim(flipdim(res.Visf,1), 2);
res.Vp = flipdim(flipdim(res.Vp,1), 2);
res.PS = flipdim(flipdim(res.PS,1), 2);
res.Ki = flipdim(flipdim(res.Ki,1), 2);

end
