
function [ori, filtered] = PerformGadgetronRecon_Perf_Debug_SpatialTemporalFiltering(debug_dir, thres_temporal, thres_spatial, wave_name, level)
% PerformGadgetronRecon_Perf_Debug_SpatialTemporalFiltering(debug_dir, thres_temporal, thres_spatial, wave_name, level)

SLC = 1;
while(isFileExist(fullfile(debug_dir, ['input_for_filter_' num2str(SLC-1) '_MAG.hdr'])))
    SLC = SLC + 1;
end

suffix = ['temporal_' num2str(thres_temporal) '_spatial_' num2str(thres_spatial) '_wave_' wave_name '_level_' num2str(level)];
ind = find(suffix=='.');
suffix(ind) = 'p';

ori = [];
filtered = [];

for slc=1:SLC-1
    data = analyze75read( fullfile(debug_dir, ['input_for_filter_' num2str(slc-1) '_MAG.hdr']) );
    res = PerformSpatioTemporalFiltering(complex(double(data)), thres_temporal, thres_spatial, wave_name, level);        
    header = CreateGtImageHeader(data, [1 1 1 1]);
    Matlab_gt_write_analyze(single(abs(res)), header, fullfile(debug_dir, ['perf_moco_' num2str(slc-1) '_' suffix '.hdr']));
    
    ori(:,:,:,slc) = data;
    filtered(:,:,:,slc) = res;
end    
