
function res = perf_load_Bulls_eye_debug_data(debugDir, save_npy, plot_flag)
% read in bulls eye data
% res = perf_load_Bulls_eye_debug_data(debugDir, save_npy)

if (nargin < 2)
    save_npy = 0;
    plot_flag = 0;
end
if (nargin < 3)
    plot_flag = 0;
end

res = struct('fmap', [], 'PS', [], 'Vb', [], 'Visf', [], 'Tc', [], 'delaymap', [], 'data', [], 'flow_data', [], ... 
    'endo_mask', [], 'epi_mask', [], 'rv_mask', [], 'rvi_mask', [], 'sectors', [], 'sectors_32', [], 'sectors_contours', [], ...
    'sectors_threshold', [], 'sectors_threshold_32', [], 'perf_for_seg', []);

res.fmap = analyze75read(fullfile(debugDir, 'Bullseye_fmap'));
res.PS = analyze75read(fullfile(debugDir, 'Bullseye_ps_map'));
res.Vb = analyze75read(fullfile(debugDir, 'Bullseye_vb_map'));
res.Visf = analyze75read(fullfile(debugDir, 'Bullseye_visf_map'));
res.Tc = analyze75read(fullfile(debugDir, 'Bullseye_tc_map'));
res.delaymap = analyze75read(fullfile(debugDir, 'Bullseye_delay_map'));

try
    res.endo_mask = analyze75read(fullfile(debugDir, 'Bullseye_endo_mask'));
    res.epi_mask = analyze75read(fullfile(debugDir, 'Bullseye_epi_mask'));
    res.rv_mask = analyze75read(fullfile(debugDir, 'Bullseye_rv_mask'));
    res.rvi_mask = analyze75read(fullfile(debugDir, 'Bullseye_rvi_mask'));
    
    res.endo_epi_rv_rvi_mask = res.epi_mask;
    res.endo_epi_rv_rvi_mask(find(res.rv_mask(:)==1)) = 3;
    res.endo_epi_rv_rvi_mask(find(res.rvi_mask(:)==1)) = 4;
    res.endo_epi_rv_rvi_mask(find(res.epi_mask(:)==1)) = 1;
    res.endo_epi_rv_rvi_mask(find(res.endo_mask(:)==1)) = 2;

catch
end

try
    res.sectors = analyze75read(fullfile(debugDir, 'Bullseye_sectors'));
    res.sectors_32 = analyze75read(fullfile(debugDir, 'Bullseye_sectors_32'));
    res.sectors_contours = analyze75read(fullfile(debugDir, 'Bullseye_sectors_contours'));
    res.sectors_threshold = analyze75read(fullfile(debugDir, 'Bullseye_sectors_threshold'));
    res.sectors_threshold_32 = analyze75read(fullfile(debugDir, 'Bullseye_sectors_threshold_32'));
catch
end

try
    res.data = analyze75read(fullfile(debugDir, 'Bullseye_data'));
    res.flow_data = analyze75read(fullfile(debugDir, 'Bullseye_flow_data'));
catch
end

res.perf_for_seg = [];
for slc=0:8
    try
        perf_seg = analyze75read(fullfile(debugDir, ['perf_seg_sax_segmentation_perf_for_seg_' num2str(slc)]));
        res.perf_for_seg = cat(4, res.perf_for_seg, perf_seg);
    catch
    end
end

if(save_npy)
    writeNPY(res.fmap, fullfile(debugDir, 'Bullseye_fmap.npy'));
    writeNPY(res.PS, fullfile(debugDir, 'Bullseye_ps_map.npy'));
    writeNPY(res.Vb, fullfile(debugDir, 'Bullseye_vb_map.npy'));
    writeNPY(res.Visf, fullfile(debugDir, 'Bullseye_visf_map.npy'));
    writeNPY(res.Tc, fullfile(debugDir, 'Bullseye_tc_map.npy'));
    writeNPY(res.delaymap, fullfile(debugDir, 'Bullseye_delay_map.npy'));
    
    if (~isempty(res.data)) writeNPY(res.data, fullfile(debugDir, 'Bullseye_data.npy')); end
    if (~isempty(res.flow_data)) writeNPY(res.flow_data, fullfile(debugDir, 'Bullseye_flow_data.npy')); end
    
    if (~isempty(res.endo_mask)) writeNPY(res.endo_mask, fullfile(debugDir, 'Bullseye_endo_mask.npy')); end
    if (~isempty(res.epi_mask)) writeNPY(res.epi_mask, fullfile(debugDir, 'Bullseye_epi_mask.npy')); end
    if (~isempty(res.rv_mask)) writeNPY(res.rv_mask, fullfile(debugDir, 'Bullseye_rv_mask.npy')); end
    if (~isempty(res.rvi_mask)) writeNPY(res.rvi_mask, fullfile(debugDir, 'Bullseye_rvi_mask.npy')); end

    if (~isempty(res.sectors)) writeNPY(res.sectors, fullfile(debugDir, 'Bullseye_sectors.npy')); end
    if (~isempty(res.sectors_32)) writeNPY(res.sectors_32, fullfile(debugDir, 'Bullseye_sectors_32.npy')); end
    if (~isempty(res.perf_for_seg)) writeNPY(res.perf_for_seg, fullfile(debugDir, 'perf_for_seg.npy')); end
end

if(plot_flag)
    figure; imagescn(res.fmap, [0 8], [1 size(res.fmap,3)], 12); PerfColorMap;
    figure; imagescn(res.perf_for_seg, [0 1.5], [1 size(res.perf_for_seg, 4)], 12, 3);
    
    if (~isempty(res.endo_mask)) figure; imagescn(res.endo_mask, [], [1 size(res.endo_mask,3)], 12); end
    if (~isempty(res.epi_mask)) figure; imagescn(res.epi_mask, [], [1 size(res.epi_mask,3)], 12); end
    if (~isempty(res.rv_mask)) figure; imagescn(res.rv_mask, [], [1 size(res.rv_mask,3)], 12); end
    if (~isempty(res.rvi_mask)) figure; imagescn(res.rvi_mask, [], [1 size(res.rvi_mask,3)], 12); end
end
