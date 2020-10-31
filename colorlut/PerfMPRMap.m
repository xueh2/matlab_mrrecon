function perf_mpr_map = PerfMPRMap(use_cmap)

if(nargin<1)
    use_cmap = 1;
end

load color_perf_mpr_map
perf_mpr_map = p;
if (use_cmap) 
    colormap(perf_mpr_map)
end