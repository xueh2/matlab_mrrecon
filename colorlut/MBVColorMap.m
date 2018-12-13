function mbv_color_map = MBVColorMap(use_cmap)

if(nargin<1)
    use_cmap = 1;
end

a = load('mbv_colormap.mat');
mbv_color_map = a.mbv_colormap;
if (use_cmap)
    colormap(mbv_color_map)
end