function [fmaps] = load_maps_for_perf(debug_dir)
% [fmaps] = load_maps_for_perf(debug_dir)

fmaps = [];

for slc=1:12
    fmap_name = fullfile(debug_dir, ['flow_map_after_hole_filling_' num2str(slc-1) '_0.hdr']);
    if exist(fmap_name)
        m = analyze75read(fmap_name);
        if slc>1 && size(m, 1) ~= size(fmaps, 1)
            m = m';
        end
        fmaps(:,:,slc) = m;
    end
end

size(fmaps)
