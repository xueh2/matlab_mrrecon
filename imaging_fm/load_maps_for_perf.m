function [fmaps] = load_maps_for_perf(debug_dir)
% [fmaps] = load_maps_for_perf(debug_dir)

fmaps = [];

for slc=1:12
    fmap_name = fullfile(debug_dir, ['flow_maps_' num2str(slc-1) '.hdr']);
    if exist(fmap_name)
        fmaps(:,:,slc) = analyze75read(fmap_name);
    end
end

size(fmaps)
