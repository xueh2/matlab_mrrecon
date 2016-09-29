function [rest, KiMap, FlowMap, aif_fig, aif, ori, moco, moco_norm, aif_ori, aif_moco] = ReviewPerfFlowResults(baseDir, resDir)
% [rest, KiMap, FlowMap, aif_fig, aif, ori, moco, moco_norm, aif_ori, aif_moco] = ReviewPerfFlowResults(baseDir, resDir)

if(nargin<2)
    resDir = 'grappa_flow_res_BTEX20_4';
end

restDir = fullfile(baseDir, resDir)
try
%         rest_perf = readGTPlusExportImageSeries_Squeeze(restDir, 105);
    r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_0.hdr'));
    r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_1.hdr'));
    r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_2.hdr'));
    rest = cat(4, r1, r2, r3);
    rest = permute(rest, [1 2 4 3]);
catch
    rest = readGTPlusExportImageSeries_Squeeze(restDir, 118);
    rest = rest/1000;
end

ori = readGTPlusExportImageSeries_Squeeze(restDir, 103);
moco = readGTPlusExportImageSeries_Squeeze(restDir, 104);
moco_norm = readGTPlusExportImageSeries_Squeeze(restDir, 105);
aif_ori = readGTPlusExportImageSeries_Squeeze(restDir, 1104);
aif_moco = readGTPlusExportImageSeries_Squeeze(restDir, 101);

a = readGTPlusExportImageSeries_Squeeze(restDir, 120);
fa = readGTPlusExportImageSeries_Squeeze(restDir, 115);
fa = fa/100;
KiMap = fa;

fa_pde = readGTPlusExportImageSeries_Squeeze(restDir, 124);
fa_pde = fa_pde/100;
FlowMap = fa_pde;

visf_a = readGTPlusExportImageSeries_Squeeze(restDir, 116);
visf_a = visf_a/10;

rest = flipdim(flipdim(rest, 2), 1);
fa = flipdim(flipdim(fa, 2), 1);
fa_pde = flipdim(flipdim(fa_pde, 2), 1);
visf_a = flipdim(flipdim(visf_a, 2), 1);

try
    aif_mask = analyze75read(fullfile(restDir, 'DebugOutput', 'AifLVMask_after_Picking.hdr'));
catch
    aif_mask = [];
end

try
    aif = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
catch
    aif = [];
end   

try
    r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'PS_maps_after_hole_filling_0.hdr'));
    r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'PS_maps_after_hole_filling_1.hdr'));
    r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'PS_maps_after_hole_filling_2.hdr'));
    PS = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
    PS = flipdim(flipdim(PS, 2), 1);
catch
end

try
    r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_0.hdr'));
    r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_1.hdr'));
    r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'blood_volume_maps_after_hole_filling_2.hdr'));
    Vp = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
    Vp = flipdim(flipdim(Vp, 2), 1);
catch
end

try
    r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'interstitial_volume_maps_0.hdr'));
    r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'interstitial_volume_maps_1.hdr'));
    r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'interstitial_volume_maps_2.hdr'));
    Visf = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
    Visf = flipdim(flipdim(Visf, 2), 1);
catch
end

try
    r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'E_maps_0.hdr'));
    r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'E_maps_1.hdr'));
    r3 = analyze75read(fullfile(restDir, 'DebugOutput', 'E_maps_2.hdr'));
    E = cat(3, r1(:,:,end), r2(:,:,end), r3(:,:,end));
    E = flipdim(flipdim(E, 2), 1);
catch
end
    
a = permute(a, [2 1 3]);
aif_fig = a;

h = figure; 
imagescn(a, [], [], 25)
saveas(h, fullfile(baseDir, [resDir '_AIF_FIG']), 'fig')

slc = size(fa, 3);
m = size(fa, 4);

figure; imagescn(fa, [0 6], [m slc], 10); PerfColorMap;

scalingFactor = 10;

h = figure('Name','Ki Maps - TwoCompExp','NumberTitle','off'); imagescn(fa(:,:,:,3), [0 6], [1 slc], scalingFactor); PerfColorMap;
saveas(h, fullfile(baseDir, [resDir '_Ki_TwoComp']), 'fig')

h = figure('Name','Gd images','NumberTitle','off'); imagescn(rest, [0 3], [], scalingFactor, 4);
saveas(h, fullfile(baseDir, [resDir '_Gd']), 'fig')

h = figure('Name','Flow Maps','NumberTitle','off'); imagescn(fa_pde(:,:,:,end), [0 6], [1 slc], scalingFactor); PerfColorMap;
saveas(h, fullfile(baseDir, [resDir '_PDE_FlowMap']), 'fig')

h = figure('Name','PDE Visf','NumberTitle','off'); imagescn(Visf(:,:,:,end), [0 100], [1 slc], scalingFactor); PerfColorMap;
saveas(h, fullfile(baseDir, [resDir '_PDE_Visf']), 'fig')

h = figure('Name','PDE PS','NumberTitle','off'); imagescn(PS(:,:,:,end), [0 10], [1 slc], scalingFactor); PerfColorMap;
saveas(h, fullfile(baseDir, [resDir '_PDE_PS']), 'fig')

h = figure('Name','PDE E','NumberTitle','off'); imagescn(E(:,:,:,end), [0 2], [1 slc], scalingFactor); PerfColorMap;
saveas(h, fullfile(baseDir, [resDir '_PDE_E']), 'fig')

h = figure('Name','PDE Vp','NumberTitle','off'); imagescn(Vp(:,:,:,end), [0 40], [1 slc], scalingFactor); PerfColorMap;
saveas(h, fullfile(baseDir, [resDir '_PDE_Vp']), 'fig')

h = figure('Name','Ori','NumberTitle','off'); imagescn(ori, [], [], scalingFactor, 4);
saveas(h, fullfile(baseDir, [resDir '_Ori']), 'fig')

h = figure('Name','MOCO','NumberTitle','off'); imagescn(moco, [], [], scalingFactor, 4);
saveas(h, fullfile(baseDir, [resDir '_MOCO']), 'fig')

h = figure('Name','MOCO NORM','NumberTitle','off'); imagescn(moco_norm, [], [], scalingFactor, 4);
saveas(h, fullfile(baseDir, [resDir '_MOCO_NORM']), 'fig')

h = figure('Name','AIF Ori','NumberTitle','off'); imagescn(aif_ori, [], [], scalingFactor, 4);
saveas(h, fullfile(baseDir, [resDir '_AIF_ORI']), 'fig')

h = figure('Name','AIF MOCO','NumberTitle','off'); imagescn(aif_moco, [], [], scalingFactor, 4);
saveas(h, fullfile(baseDir, [resDir '_AIF_MOCO']), 'fig')

if(~isempty(aif_mask))
    N = size(aif_moco, 4);
    
    ind = find(aif_mask>0);
    
    for n=1:N
        a2D = aif_moco(:,:,1, n);
        a2D(ind(:)) = 1000;
        aif_moco(:,:,1, n) = a2D;
    end
    
    h = figure('Name','AIF MOCO with mask','NumberTitle','off'); imagescn(aif_moco, [], [], scalingFactor, 4);
    saveas(h, fullfile(baseDir, [resDir '_AIF_MOCO_With_Mask']), 'fig')
end

delta = 0.5;

if(~isempty(aif))
    h = figure
    hold on
    plot(delta*[0:numel(aif)-1], aif, 'LineWidth',2)
    hold off
    set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
    box on
    grid on
    grid MINOR
    xlabel('second')
    ylabel('Gd [mmol/ml]')
    title('AIF in Gd')

    saveas(h, fullfile(baseDir, [resDir '_AIF_Curves']), 'fig')
end
