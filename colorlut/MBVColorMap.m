function mbv_color_map = MBVColorMap

a = load('mbv_colormap.mat');
mbv_color_map = a.mbv_colormap;
colormap(mbv_color_map)