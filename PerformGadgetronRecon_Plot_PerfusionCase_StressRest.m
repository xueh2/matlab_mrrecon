
function PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir, stressCase, restCase, onlyReview, baseDir)
% PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir, stressCase, onlyReview, restCase)
% PerformGadgetronRecon_Plot_PerfusionCase_StressRest('I:\ReconResults\BARTS', stressCase, restCase, onlyReview)

if(nargin < 4)
    onlyReview = 0;
end

if(nargin < 5)
    baseDir = resDir;
end

[configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(stressCase);

if(~onlyReview)   
    try
        restDir = fullfile(resDir, study_dates, restCase)
        try
    %         rest_perf = readGTPlusExportImageSeries_Squeeze(restDir, 105);
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_2.hdr'));
            rest_perf = cat(4, r1, r2, r3);
            rest_perf = permute(rest_perf, [1 2 4 3]);
            rest_perf = flipdim(rest_perf, 2);
        catch
            rest_perf = readGTPlusExportImageSeries_Squeeze(restDir, 111);
        end

        r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'perf_0.hdr'));
        r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'perf_1.hdr'));
        r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'perf_2.hdr'));
        ori_rest = flipdim(cat(4, r1, r2, r3), 2);
        
        r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'moco_0_MAG.hdr'));
        r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'moco_1_MAG.hdr'));
        r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'moco_2_MAG.hdr'));
        moco_rest = flipdim(cat(4, r1, r2, r3), 2);

        r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'SRNorm_0.hdr'));
        r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'SRNorm_1.hdr'));
        r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'SRNorm_2.hdr'));
        moco_norm_rest = flipdim(cat(4, r1, r2, r3), 2);

        r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_moco.hdr'));
        r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_moco_second_echo.hdr'));
        aif_moco_rest = flipdim(cat(4, r1, r2), 2);

        try
            aif_rest = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
        catch
            aif_rest = [];
        end

        try
            aif_rest_baseline_corrected = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_cin_echo0_all_signal_baseline_corrected.hdr'));
        catch
            aif_rest_baseline_corrected = [];
        end
        
        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'flow_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'flow_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'flow_maps_after_hole_filling_2.hdr'));
            flow = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            flow_rest = flipdim(flow, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'Ki_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'Ki_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'Ki_maps_after_hole_filling_2.hdr'));
            Ki = cat(4, r1, r2, r3);
            Ki_rest = flipdim(Ki, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'PS_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'PS_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'PS_maps_after_hole_filling_2.hdr'));
            PS = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            PS_rest = flipdim(PS, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_2.hdr'));
            Vp = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            Vp_rest = flipdim(Vp, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'interstitial_volume_maps_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'interstitial_volume_maps_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'interstitial_volume_maps_2.hdr'));
            Visf = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            Visf_rest = flipdim(Visf, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'E_maps_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'E_maps_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'E_maps_2.hdr'));
            E = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            E_rest = flipdim(E, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'BTEX_SD_maps_0_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'BTEX_SD_maps_1_0.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'BTEX_SD_maps_2_0.hdr'));
            SDMap = cat(3, r1, r2, r3);
            SDMap_rest = flipdim(SDMap, 2);
        catch
        end
        
        has_rest = 1;
       
    catch
        has_rest = 0;
        rest_perf = 0;
        fa = 0;
        a = 0;
        aif_rest = 0;
    end

    try
        stressDir = fullfile(resDir, study_dates, stressCase)
        try
            % stress_perf = readGTPlusExportImageSeries_Squeeze(stressDir, 105);
            r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_0.hdr'));
            r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_1.hdr'));
            r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_2.hdr'));
            stress_perf = cat(4, r1, r2, r3);
            stress_perf = permute(stress_perf, [1 2 4 3]);
            stress_perf = flipdim(stress_perf, 2);
        catch
            stress_perf = readGTPlusExportImageSeries_Squeeze(stressDir, 111);
        end

        r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'perf_0.hdr'));
        r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'perf_1.hdr'));
        r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'perf_2.hdr'));
        ori_stress = cat(4, r1, r2, r3);
        ori_stress = flipdim(ori_stress, 2);

        r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'moco_0_MAG.hdr'));
        r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'moco_1_MAG.hdr'));
        r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'moco_2_MAG.hdr'));
        moco_stress = cat(4, r1, r2, r3);
        moco_stress = flipdim(moco_stress, 2);

        r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'SRNorm_0.hdr'));
        r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'SRNorm_1.hdr'));
        r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'SRNorm_2.hdr'));
        moco_norm_stress = cat(4, r1, r2, r3);
        moco_norm_stress = flipdim(moco_norm_stress, 2);
        
        r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'aif_moco.hdr'));
        r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'aif_moco_second_echo.hdr'));
        aif_moco_stress = cat(4, r1, r2);
        aif_moco_stress = flipdim(aif_moco_stress, 2);
        
        try
            aif_stress = analyze75read(fullfile(stressDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
        catch
            aif_stress = [];
        end

        try
            aif_stress_baseline_corrected = analyze75read(fullfile(stressDir, 'DebugOutput', 'aif_cin_echo0_all_signal_baseline_corrected.hdr'));
        catch
            aif_stress_baseline_corrected = [];
        end
        
        try
            r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'flow_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'flow_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'flow_maps_after_hole_filling_2.hdr'));
            flow = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            flow_stress = flipdim(flow, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'Ki_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'Ki_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'Ki_maps_after_hole_filling_2.hdr'));
            Ki = cat(4, r1, r2, r3);
            Ki_stress = flipdim(Ki, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'PS_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'PS_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'PS_maps_after_hole_filling_2.hdr'));
            PS = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            PS_stress = flipdim(PS, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_2.hdr'));
            Vp = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            Vp_stress = flipdim(Vp, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'interstitial_volume_maps_0.hdr'));
            r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'interstitial_volume_maps_1.hdr'));
            r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'interstitial_volume_maps_2.hdr'));
            Visf = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            Visf_stress = flipdim(Visf, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'E_maps_0.hdr'));
            r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'E_maps_1.hdr'));
            r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'E_maps_2.hdr'));
            E = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            E_stress = flipdim(E, 2);
        catch
        end

        try
            r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'BTEX_SD_maps_0_0.hdr'));
            r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'BTEX_SD_maps_1_0.hdr'));
            r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'BTEX_SD_maps_2_0.hdr'));
            SDMap = cat(3, r1, r2, r3);
            SDMap_stress = flipdim(SDMap, 2);
        catch
        end
        
        has_stress = 1;
    catch
        has_stress = 0;
        stress_perf = 0;
        fb = 0;
        b = 0;
        aif_stress = 0;
    end
end

% scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time
resDir = [scannerID '_' patientID '_' studyID '_' study_dates];

figDir = fullfile(baseDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' resDir '_Figure']);
disp(figDir)

if(~onlyReview)
    if(has_stress || has_rest)  
               
        if(exist(figDir)==7)
            try
                rmdir(figDir, 's');
            catch
            end
        end
        mkdir(figDir);
        
        save(fullfile(figDir, 'rest.mat'), 'rest_perf', 'ori_rest', 'moco_rest', 'moco_norm_rest', 'aif_moco_rest', 'aif_rest', 'aif_rest_baseline_corrected', 'flow_rest', 'Ki_rest', 'PS_rest', 'Vp_rest', 'Visf_rest', 'E_rest', 'SDMap_rest');
        save(fullfile(figDir, 'stress.mat'), 'stress_perf', 'ori_stress', 'moco_stress', 'moco_norm_stress', 'aif_moco_stress', 'aif_stress', 'aif_stress_baseline_corrected', 'flow_stress', 'Ki_stress', 'PS_stress', 'Vp_stress', 'Visf_stress', 'E_stress', 'SDMap_stress');
    end
else
    if(exist(figDir)~=7)
        disp([figDir ' - does not exist']);
        return;
    end
    
    [nameFigs, numFigs] = findFILE(figDir, '*.fig');
    
    if(numFigs > 0 )
        has_stress = 1;
        has_rest = 1;
    else
        disp([figDir ' - no figures are found']);
        return;
    end
end

if(has_stress & has_rest)
%     a = permute(a, [2 1 3]);
%     b = permute(b, [2 1 3]);
% 
%     h = figure; 
%     imagescn(cat(4, a, b), [], [], 25)
%     saveas(h, fullfile(figDir, [resDir '_AIF_FIG']), 'fig')

    if(~onlyReview)
        slc = size(Ki_rest, 4);
        m = size(Ki_rest, 3);

        if(size(Ki_rest, 1)~=size(Ki_stress,1) || size(Ki_rest, 2)~=size(Ki_stress,2) )
            disp(['Image size mismatch - Rest :' num2str(size(Ki_rest)) ' - Stress : ' num2str(size(Ki_stress))]);
            return;
        end

        scalingFactor = 10;

        figure; imagescn(cat(4, Ki_stress, Ki_rest), [0 6], [m 2*slc], scalingFactor); PerfColorMap;
    else
        slc = 3;
    end
    scrsz = get(0, 'ScreenSize');
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_Ki_TwoComp']);
    if(onlyReview)
        openfig(figName);
        
        for s=1:slc
            figName = fullfile(figDir, [resDir '_PDE_StressKiMap - ' num2str(s)]);
            openfig(figName);
        end
        
        for s=1:slc
            figName = fullfile(figDir, [resDir '_PDE_RestKiMap - ' num2str(s)]);
            openfig(figName);
        end
        
    else
        h = figure('Name','Ki Maps - TwoCompExp','NumberTitle','off'); imagescn(cat(4, Ki_stress(:,:,3,:), Ki_rest(:,:,3,:)), [0 6], [2 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
        
        for s=1:slc
            h = figure('Name',['Stress Ki maps - ' num2str(s)],'NumberTitle','off'); imagescn(Ki_stress(:,:,3,s), [0 6], [1 1], scalingFactor); PerfColorMap;
            figName = fullfile(figDir, [resDir '_PDE_StressKiMap - ' num2str(s)]);
            saveas(h, figName, 'fig')
        end
        
        for s=1:slc
            h = figure('Name',['Rest Ki maps - ' num2str(s)],'NumberTitle','off'); imagescn(Ki_rest(:,:,3,s), [0 6], [1 1], scalingFactor); PerfColorMap;
            figName = fullfile(figDir, [resDir '_PDE_RestKiMap - ' num2str(s)]);
            saveas(h, figName, 'fig')
        end
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_Gd']);
    if(onlyReview)
        % openfig(figName);
    else
        NS = size(stress_perf, 4);
        NR = size(rest_perf, 4);
        
        NN = NS;
        if(NN>NR) NN = NR; end
        
        h = figure('Name','Rest-Stress Gd','NumberTitle','off'); imagescn(cat(3, stress_perf(:,:,:,1:NN), rest_perf(:,:,:,1:NN)), [0 3], [2 slc], 10, 4);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_FlowMap']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Flow maps','NumberTitle','off'); imagescn(cat(3, flow_stress(:,:,:,end), flow_rest(:,:,:,end)), [0 6], [2 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_PDE_StressFlowMap']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Stress Flow maps','NumberTitle','off'); imagescn(flow_stress(:,:,:,end), [0 6], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
        
        for s=1:slc
            h = figure('Name',['Stress Flow maps - ' num2str(s)],'NumberTitle','off'); imagescn(flow_stress(:,:,s,end), [0 6], [1 1], scalingFactor); PerfColorMap;
            figName = fullfile(figDir, [resDir '_PDE_StressFlowMap - ' num2str(s)]);
            saveas(h, figName, 'fig')
        end
    end
    
    figName = fullfile(figDir, [resDir '_PDE_RestFlowMap']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Rest Flow maps','NumberTitle','off'); imagescn(flow_rest(:,:,:,end), [0 6], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
        
        for s=1:slc
            h = figure('Name',['Rest Flow maps - ' num2str(s)],'NumberTitle','off'); imagescn(flow_rest(:,:,s,end), [0 6], [1 1], scalingFactor); PerfColorMap;
            figName = fullfile(figDir, [resDir '_PDE_RestFlowMap - ' num2str(s)]);
            saveas(h, figName, 'fig')
        end
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_Visf']);
    if(onlyReview)
        % openfig(figName);
    else
        h = figure('Name','PDE Visf','NumberTitle','off'); imagescn(cat(3, Visf_stress(:,:,:,end), Visf_rest(:,:,:,end)), [0 100], [2 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_PS']);
    if(onlyReview)
        % openfig(figName);
    else
        h = figure('Name','PDE PS','NumberTitle','off'); imagescn(cat(3, PS_stress(:,:,:,end), PS_rest(:,:,:,end)), [0 10], [2 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_E']);
    if(onlyReview)
        % openfig(figName);
    else
        h = figure('Name','PDE E','NumberTitle','off');; imagescn(cat(3, E_stress(:,:,:,end), E_rest(:,:,:,end)), [0 2], [2 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_Vp']);
    if(onlyReview)
        % openfig(figName);
    else
        h = figure('Name','PDE Vp','NumberTitle','off');; imagescn(cat(3, Vp_stress(:,:,:,end), Vp_rest(:,:,:,end)), [0 20], [2 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_Flow_SD']);
    if(onlyReview)
        % openfig(figName);
    else
        h = figure('Name','PDE Flow SD','NumberTitle','off');; imagescn(cat(3, SDMap_stress(:,:,:,end), SDMap_rest(:,:,:,end)), [0 0.8], [2 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_Ori']);
    if(onlyReview)
        % openfig(figName);
    else
        h = figure('Name','Rest-Stress Original','NumberTitle','off'); imagescn(cat(4, ori_stress(:,:,1:NN,:), ori_rest(:,:,1:NN,:)), [], [2 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_MOCO']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Rest-Stress MOCO','NumberTitle','off'); imagescn(cat(4, moco_stress(:,:,1:NN,:), moco_rest(:,:,1:NN,:)), [], [2 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_CoMOCO']);
    if(onlyReview)
        openfig(figName);
    else
        ss = moco_stress(:,:,[1 4 NN],:);
        rr = moco_rest(:,:,[1 4 NN],:);
        
        h = figure('Name','Rest-Stress SR-PD CoMOCO','NumberTitle','off'); imagescn(cat(4, ss, rr), [], [2*slc 3], scalingFactor);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_MOCO_NORM']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Rest-Stress MOCO NORM','NumberTitle','off'); imagescn(cat(4, moco_norm_stress(:,:,1:NN,:), moco_norm_rest(:,:,1:NN,:)), [0 3], [2 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Stress_Rest_AIF_MOCO']);
    if(onlyReview)
        % openfig(figName);
    else
        h = figure('Name','AIF MOCO','NumberTitle','off'); imagescn(cat(4, aif_moco_stress(:,:,1:NN,:), aif_moco_rest(:,:,1:NN,:)), [], [], 10, 3);
        saveas(h, figName, 'fig')
    end
    
    delta = 0.5;

    figName = fullfile(figDir, [resDir '_AIF_Stress_Rest_Curves']);

    if(onlyReview)
        openfig(figName);
    else        
        h = figure
        hold on
        plot(delta*[0:numel(aif_stress)-1], aif_stress, 'r', 'LineWidth',2);
        plot(delta*[0:numel(aif_rest)-1], aif_rest, 'LineWidth',2)
        hold off
        set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
        box on
        grid on
        grid MINOR
        legend('stress', 'rest')
        xlabel('second')
        ylabel('Gd [mmol/ml]')
        title('AIF in Gd')

        saveas(h, figName, 'fig')
    end
else
    if(has_stress & ~has_rest)
%         b = permute(b, [2 1 3]);
% 
%         h = figure; 
%         imagescn(b, [], [], 25)
%         saveas(h, fullfile(figDir, [resDir '_AIF_Stress_FIG']), 'fig')

        slc = size(fb, 3);
        m = size(fb, 4);

        figure; imagescn(fb, [0 6], [m slc], scalingFactor); PerfColorMap;

        figure; imagescn(stress_perf, [0 3], [], 10, 4);

        delta = 0.5;

        if(~isempty(aif_stress))
            figure
            hold on
            plot(delta*[0:numel(aif_stress)-1], aif_stress, 'r', 'LineWidth',2);
            hold off
            set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
            box on
            grid on
            grid MINOR
            legend('stress')
            xlabel('second')
            ylabel('Gd [mmol/ml]')
            title('AIF in Gd')
        end
    else
        if(~has_stress & has_rest)
%             a = permute(a, [2 1 3]);
% 
%             h = figure; 
%             imagescn(a, [], [], 25)
%             saveas(h, fullfile(figDir, [resDir '_AIF_Rest_FIG']), 'fig')

            slc = size(fa, 3);
            m = size(fa, 4);

            figure; imagescn(fa, [0 6], [m slc], scalingFactor); PerfColorMap;

            figure; imagescn(rest_perf, [0 3], [], 10, 4);

            delta = 0.5;

            if(~isempty(aif_rest))
                figure
                hold on
                plot(delta*[0:numel(aif_rest)-1], aif_rest, 'r', 'LineWidth',2);
                hold off
                set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
                box on
                grid on
                grid MINOR
                legend('rest')
                xlabel('second')
                ylabel('Gd [mmol/ml]')
                title('AIF in Gd')
            end
        end
    end
end
