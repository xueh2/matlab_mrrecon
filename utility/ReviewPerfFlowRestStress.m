function [rest_perf, fa, fa_pde, a, aif_rest, stress_perf, fb, fb_pde, b, PS_rest, PS_stress, E_rest, E_stress, Vp_rest, Vp_stress, Visf_rest, Visf_stress, aif_stress, ori_rest, ori_stress, moco_rest, moco_stress, moco_norm_rest, moco_norm_stress, aif_ori_rest, aif_moco_rest, aif_ori_stress, aif_moco_stress] = ReviewPerfFlowRestStress(baseDir, resDir)
% [rest_perf, fa, a, aif_rest, stress_perf, fb, b, PS_rest, PS_stress, E_rest, E_stress, Vp_rest, Vp_stress, Visf_rest, Visf_stress, aif_stress, ori_rest, ori_stress, moco_rest, moco_stress, moco_norm_rest, moco_norm_stress, aif_ori_rest, aif_moco_rest, aif_ori_stress, aif_moco_stress] = ReviewPerfFlowRestStress(baseDir, resDir)

if(nargin<2)
    resDir = 'grappa_flow_res_new3';
end

try
    [subdirs, nums] = FindSubDirs(fullfile(baseDir, 'rest'))    
    if ( ~isempty(strfind( lower(subdirs{1}), 'mini')) | ~isempty(strfind(subdirs{1}, 'MINI')) )
        rest = subdirs{2};
    else
        rest = subdirs{1};
    end

    restDir = fullfile(baseDir, 'rest', rest, resDir)
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
    
    ori_rest = readGTPlusExportImageSeries_Squeeze(restDir, 103);
    moco_rest = readGTPlusExportImageSeries_Squeeze(restDir, 104);
    moco_norm_rest = readGTPlusExportImageSeries_Squeeze(restDir, 105);
    aif_ori_rest = readGTPlusExportImageSeries_Squeeze(restDir, 1104);
    aif_moco_rest = readGTPlusExportImageSeries_Squeeze(restDir, 101);
    
    a = readGTPlusExportImageSeries_Squeeze(restDir, 120);
    fa = readGTPlusExportImageSeries_Squeeze(restDir, 115);
    fa = fa/100;
    
    fa_pde = readGTPlusExportImageSeries_Squeeze(restDir, 124);
    fa_pde = fa_pde/100;
    
    visf_a = readGTPlusExportImageSeries_Squeeze(restDir, 116);
    visf_a = visf_a/10;
    
    rest_perf = flipdim(flipdim(rest_perf, 2), 1);
    fa = flipdim(flipdim(fa, 2), 1);
    fa_pde = flipdim(flipdim(fa_pde, 2), 1);
    visf_a = flipdim(flipdim(visf_a, 2), 1);
    
    try
        aif_rest = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
    catch
        aif_rest = [];
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

try
    [subdirs, nums] = FindSubDirs(fullfile(baseDir, 'stress'))
    if ( ~isempty(strfind( lower(subdirs{1}), 'mini')) | ~isempty(strfind(subdirs{1}, 'MINI')) )
        stress = subdirs{2};
    else
        stress = subdirs{1};
    end

    stressDir = fullfile(baseDir, 'stress', stress, resDir)
    try
        % stress_perf = readGTPlusExportImageSeries_Squeeze(stressDir, 105);
        r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_0.hdr'));
        r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_1.hdr'));
        r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_2.hdr'));
        stress_perf = cat(4, r1, r2, r3);
        stress_perf = permute(stress_perf, [1 2 4 3]);
    catch
        stress_perf = readGTPlusExportImageSeries_Squeeze(stressDir, 111);
    end
    
    ori_stress = readGTPlusExportImageSeries_Squeeze(stressDir, 103);
    moco_stress = readGTPlusExportImageSeries_Squeeze(stressDir, 104);
    moco_norm_stress = readGTPlusExportImageSeries_Squeeze(stressDir, 105);
    aif_ori_stress = readGTPlusExportImageSeries_Squeeze(stressDir, 1104);
    aif_moco_stress = readGTPlusExportImageSeries_Squeeze(stressDir, 101);
    
    b = readGTPlusExportImageSeries_Squeeze(stressDir, 120);
    fb = readGTPlusExportImageSeries_Squeeze(stressDir, 115);
    fb = fb/100;
    fb_pde = readGTPlusExportImageSeries_Squeeze(stressDir, 124);
    fb_pde = fb_pde/100;

    visf_b = readGTPlusExportImageSeries_Squeeze(stressDir, 116);
    visf_b = visf_b/10;

    stress_perf = flipdim(flipdim(stress_perf, 2), 1);
    fb = flipdim(flipdim(fb, 2), 1);
    fb_pde = flipdim(flipdim(fb_pde, 2), 1);
    visf_b = flipdim(flipdim(visf_b, 2), 1);    
    
    try
        aif_stress = analyze75read(fullfile(stressDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
    catch
        aif_stress = [];
    end
    
    try
        r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'PS_maps_after_hole_filling_0.hdr'));
        r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'PS_maps_after_hole_filling_1.hdr'));
        r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'PS_maps_after_hole_filling_2.hdr'));
        PS = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
        PS = flipdim(flipdim(PS, 2), 1);
        PS_stress = PS;
    catch
    end
    
    try
        r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_0.hdr'));
        r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_1.hdr'));
        r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_2.hdr'));
        Vp = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
        Vp = flipdim(flipdim(Vp, 2), 1);
        Vp_stress = Vp;
    catch
    end
    
    try
        r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'interstitial_volume_maps_0.hdr'));
        r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'interstitial_volume_maps_1.hdr'));
        r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'interstitial_volume_maps_2.hdr'));
        Visf = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
        Visf = flipdim(flipdim(Visf, 2), 1);
        Visf_stress = Visf;
    catch
    end
    
    try
        r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'E_maps_0.hdr'));
        r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'E_maps_1.hdr'));
        r3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'E_maps_2.hdr'));
        E = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
        E = flipdim(flipdim(E, 2), 1);
        E_stress = E;
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

if(has_stress & has_rest)
    a = permute(a, [2 1 3]);
    b = permute(b, [2 1 3]);

    h = figure; 
    imagescn(cat(4, a, b), [], [], 25)
    saveas(h, fullfile(baseDir, [resDir '_AIF_FIG']), 'fig')

    slc = size(fa, 3);
    m = size(fa, 4);

    scalingFactor = 10;
    
    figure; imagescn(cat(3, fa, fb), [0 6], [m 2*slc], scalingFactor); PerfColorMap;
    
    h = figure('Name','Ki Maps - TwoCompExp','NumberTitle','off'); imagescn(cat(3, fa(:,:,:,3), fb(:,:,:,3)), [0 6], [2 slc], scalingFactor); PerfColorMap;
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_Ki_TwoComp']), 'fig')
    
    h = figure('Name','Rest-Stress Gd','NumberTitle','off'); imagescn(cat(3, rest_perf, stress_perf), [0 3], [], 10, 4);
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_Gd']), 'fig')

    h = figure('Name','Flow maps','NumberTitle','off'); imagescn(cat(3, fa_pde(:,:,:,end), fb_pde(:,:,:,end)), [0 6], [2 slc], scalingFactor); PerfColorMap;
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_PDE_FlowMap']), 'fig')
       
    h = figure('Name','PDE Visf','NumberTitle','off'); imagescn(cat(3, Visf_rest(:,:,:,end), Visf_stress(:,:,:,end)), [0 100], [2 slc], scalingFactor); PerfColorMap;
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_PDE_Visf']), 'fig')
    
    h = figure('Name','PDE PS','NumberTitle','off'); imagescn(cat(3, PS_rest(:,:,:,end), PS_stress(:,:,:,end)), [0 10], [2 slc], scalingFactor); PerfColorMap;
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_PDE_PS']), 'fig')
    
    h = figure('Name','PDE E','NumberTitle','off');; imagescn(cat(3, E_rest(:,:,:,end), E_stress(:,:,:,end)), [0 2], [2 slc], scalingFactor); PerfColorMap;
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_PDE_E']), 'fig')
    
    h = figure('Name','PDE Vp','NumberTitle','off');; imagescn(cat(3, Vp_rest(:,:,:,end), Vp_stress(:,:,:,end)), [0 20], [2 slc], scalingFactor); PerfColorMap;
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_PDE_Vp']), 'fig')
    
    h = figure('Name','Rest-Stress Original','NumberTitle','off'); imagescn(cat(3, ori_rest, ori_stress), [], [], scalingFactor, 4);
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_Ori']), 'fig')
    
    h = figure('Name','Rest-Stress MOCO','NumberTitle','off'); imagescn(cat(3, moco_rest, moco_stress), [], [], scalingFactor, 4);
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_MOCO']), 'fig')
    
    h = figure('Name','Rest-Stress MOCO NORM','NumberTitle','off'); imagescn(cat(3, moco_norm_rest, moco_norm_stress), [], [], scalingFactor, 4);
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_MOCO_NORM']), 'fig')
    
    h = figure('Name','AIF Original','NumberTitle','off'); imagescn(cat(3, aif_ori_rest, aif_ori_stress), [], [], 10, 4);
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_AIF_ORI']), 'fig')
    
    h = figure('Name','AIF MOCO','NumberTitle','off'); imagescn(cat(3, aif_moco_rest, aif_moco_stress), [], [], 10, 4);
    saveas(h, fullfile(baseDir, [resDir '_Rest_Stress_AIF_MOCO']), 'fig')
    
    delta = 0.5;

    if(~isempty(aif_rest) & ~isempty(aif_stress) )
        h = figure
        hold on
        plot(delta*[0:numel(aif_rest)-1], aif_rest, 'LineWidth',2)
        plot(delta*[0:numel(aif_stress)-1], aif_stress, 'r', 'LineWidth',2);
        hold off
        set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
        box on
        grid on
        grid MINOR
        legend('rest', 'stress')
        xlabel('second')
        ylabel('Gd [mmol/ml]')
        title('AIF in Gd')
        
        saveas(h, fullfile(baseDir, [resDir '_AIF_Rest_Stress_Curves']), 'fig')
    end
else
    if(has_stress & ~has_rest)
        b = permute(b, [2 1 3]);

        h = figure; 
        imagescn(b, [], [], 25)
        saveas(h, fullfile(baseDir, [resDir '_AIF_FIG']), 'fig')

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
            a = permute(a, [2 1 3]);

            h = figure; 
            imagescn(a, [], [], 25)
            saveas(h, fullfile(baseDir, [resDir '_AIF_FIG']), 'fig')

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
