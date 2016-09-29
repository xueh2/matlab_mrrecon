
function PerformGadgetronRecon_Plot_PerfusionCase(resDir, restCase, onlyReview, baseDir)
% PerformGadgetronRecon_Plot_PerfusionCase(resDir, restCase, onlyReview, baseDir)
% PerformGadgetronRecon_Plot_PerfusionCase('I:\ReconResults\ROYALFREE', restCase, onlyReview)

if(nargin < 3)
    onlyReview = 0;
end

if(nargin < 4)
    baseDir = resDir;
end

[configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(restCase);

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
        catch
            rest_perf = readGTPlusExportImageSeries_Squeeze(restDir, 111);
        end

        r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'perf_0.hdr'));
        r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'perf_1.hdr'));
        r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'perf_2.hdr'));
        ori_rest = cat(4, r1, r2, r3);

        r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'moco_0_MAG.hdr'));
        r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'moco_1_MAG.hdr'));
        r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'moco_2_MAG.hdr'));
        moco_rest = cat(4, r1, r2, r3);

        r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'SRNorm_0.hdr'));
        r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'SRNorm_1.hdr'));
        r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'SRNorm_2.hdr'));
        moco_norm_rest = cat(4, r1, r2, r3);

        r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_moco.hdr'));
        r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_moco_second_echo.hdr'));
        aif_moco_rest = cat(4, r1, r2);

        %     ori_rest = readGTPlusExportImageSeries_Squeeze(restDir, 103);
    %     moco_rest = readGTPlusExportImageSeries_Squeeze(restDir, 104);
    %     moco_norm_rest = readGTPlusExportImageSeries_Squeeze(restDir, 105);
    %     aif_ori_rest = readGTPlusExportImageSeries_Squeeze(restDir, 1104);
        % aif_moco_rest = readGTPlusExportImageSeries_Squeeze(restDir, 101);

    %     a = readGTPlusExportImageSeries_Squeeze(restDir, 120);
    %     fa = readGTPlusExportImageSeries_Squeeze(restDir, 115);
    %     fa = fa/100;
    %     
    %     fa_pde = readGTPlusExportImageSeries_Squeeze(restDir, 124);
    %     fa_pde = fa_pde/100;
    %     
    %     visf_a = readGTPlusExportImageSeries_Squeeze(restDir, 116);
    %     visf_a = visf_a/10;

        rest_perf = flipdim(flipdim(rest_perf, 2), 1);
    %     fa = flipdim(flipdim(fa, 2), 1);
    %     fa_pde = flipdim(flipdim(fa_pde, 2), 1);
    %     visf_a = flipdim(flipdim(visf_a, 2), 1);

        try
            aif_rest = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
        catch
            aif_rest = [];
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'flow_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'flow_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'flow_maps_after_hole_filling_2.hdr'));
            flow = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            flow = flipdim(flipdim(flow, 2), 1);
            flow_rest = flow;
        catch
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'Ki_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'Ki_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'Ki_maps_after_hole_filling_2.hdr'));
            Ki = cat(4, r1, r2, r3);
            Ki = flipdim(flipdim(Ki, 2), 1);
            Ki_rest = Ki;
        catch
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'PS_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'PS_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'PS_maps_after_hole_filling_2.hdr'));
            PS = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            PS = flipdim(flipdim(PS, 2), 1);
            PS_rest = PS;
        catch
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_2.hdr'));
            Vp = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            Vp = flipdim(flipdim(Vp, 2), 1);
            Vp_rest = Vp;
        catch
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'interstitial_volume_maps_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'interstitial_volume_maps_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'interstitial_volume_maps_2.hdr'));
            Visf = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            Visf = flipdim(flipdim(Visf, 2), 1);
            Visf_rest = Visf;
        catch
        end

        try
            r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'E_maps_0.hdr'));
            r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'E_maps_1.hdr'));
            r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'E_maps_2.hdr'));
            E = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
            E = flipdim(flipdim(E, 2), 1);
            E_rest = E;
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
end

% scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time
resDir = [scannerID '_' patientID '_' studyID '_' study_dates];

figDir = fullfile(baseDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' resDir '_Figure']);

if(~onlyReview)
    if(has_rest)        
        if(exist(figDir)==7)
            rmdir(figDir, 's');
        end
        mkdir(figDir);
    end
else
    has_rest = 1;
end

if(has_rest)
    if(~onlyReview)
        slc = size(Ki_rest, 4);
        m = size(Ki_rest, 3);

        scalingFactor = 10;

        figure; imagescn(Ki_rest, [0 6], [m slc], scalingFactor); PerfColorMap;
    end
    scrsz = get(0, 'ScreenSize');
    
    figName = fullfile(figDir, [resDir '_Rest_Ki_TwoComp']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Ki Maps - TwoCompExp','NumberTitle','off'); imagescn(Ki_rest(:,:,3,:), [0 6], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_Gd']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Rest Gd','NumberTitle','off'); imagescn(rest_perf, [0 3], [1 slc], 10, 4);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_FlowMap']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Flow maps','NumberTitle','off'); imagescn(flow_rest(:,:,:,end), [0 6], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_Visf']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','PDE Visf','NumberTitle','off'); imagescn(Visf_rest(:,:,:,end), [0 100], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_PS']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','PDE PS','NumberTitle','off'); imagescn(PS_rest(:,:,:,end), [0 10], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_E']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','PDE E','NumberTitle','off');; imagescn(E_rest(:,:,:,end), [0 2], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_Vp']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','PDE Vp','NumberTitle','off');; imagescn(Vp_rest(:,:,:,end), [0 20], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_Ori']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Rest Original','NumberTitle','off'); imagescn(ori_rest, [], [1 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_MOCO']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Rest MOCO','NumberTitle','off'); imagescn(moco_rest, [], [1 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_MOCO_NORM']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Rest MOCO NORM','NumberTitle','off'); imagescn(moco_norm_rest, [], [1 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_AIF_MOCO']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','AIF MOCO','NumberTitle','off'); imagescn(aif_moco_rest, [], [], 10, 3);
        saveas(h, figName, 'fig')
    end
    
    delta = 0.5;

    figName = fullfile(figDir, [resDir '_AIF_Rest_Curves']);

    if(onlyReview)
        openfig(figName);
    else        
        h = figure
        hold on
        plot(delta*[0:numel(aif_rest)-1], aif_rest, 'LineWidth',2)
        hold off
        set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
        box on
        grid on
        grid MINOR
        legend('rest')
        xlabel('second')
        ylabel('Gd [mmol/ml]')
        title('AIF in Gd')

        saveas(h, figName, 'fig')
    end
end
