function [rest, KiMap, FlowMap, aif_fig, aif, ori, moco, moco_norm, aif_ori, aif_moco] = ReviewPerfFlowNonLinearResults(baseDir, resDir)
% [rest, KiMap, FlowMap, aif_fig, aif, ori, moco, moco_norm, aif_ori, aif_moco] = ReviewPerfFlowNonLinearResults(baseDir, resDir)

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

ori = readGTPlusExportImageSeries_Squeeze(restDir, 109);
moco = readGTPlusExportImageSeries_Squeeze(restDir, 109);
moco_norm = readGTPlusExportImageSeries_Squeeze(restDir, 111);
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
    aif = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
catch
    aif = [];
end   

a = permute(a, [2 1 3]);
aif_fig = a;

h = figure; 
imagescn(a, [], [], 25)
saveas(h, fullfile(baseDir, [resDir '_AIF_FIG']), 'fig')

slc = size(fa, 3);
m = size(fa, 4);

figure; imagescn(fa, [0 6], [m slc], 10); PerfColorMap;

h = figure; imagescn(fa(:,:,:,3), [0 6], [1 slc], 10); PerfColorMap;
saveas(h, fullfile(baseDir, [resDir '_Ki_TwoComp']), 'fig')

h = figure; imagescn(rest, [0 3], [], 10, 4);
saveas(h, fullfile(baseDir, [resDir '_Gd']), 'fig')

h = figure; imagescn(fa_pde(:,:,:,end), [0 6], [1 slc], 10); PerfColorMap;
saveas(h, fullfile(baseDir, [resDir '_FlowMap_TwoComp']), 'fig')

h = figure; imagescn(visf_a(:,:,:,end), [0 100], [1 slc], 10); PerfColorMap;
saveas(h, fullfile(baseDir, [resDir '_PDE_Visf']), 'fig')

h = figure; imagescn(ori, [], [], 10, 4);
saveas(h, fullfile(baseDir, [resDir '_Ori']), 'fig')
h = figure; imagescn(moco, [], [], 10, 4);
saveas(h, fullfile(baseDir, [resDir '_MOCO']), 'fig')
h = figure; imagescn(moco_norm, [], [], 10, 4);
saveas(h, fullfile(baseDir, [resDir '_MOCO_NORM']), 'fig')

h = figure; imagescn(aif_ori, [], [], 10, 4);
saveas(h, fullfile(baseDir, [resDir '_AIF_ORI']), 'fig')

h = figure; imagescn(aif_moco, [], [], 10, 4);
saveas(h, fullfile(baseDir, [resDir '_AIF_MOCO']), 'fig')

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
