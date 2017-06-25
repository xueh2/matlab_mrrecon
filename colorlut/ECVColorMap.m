function ecv_color_map = ECVColorMap

a = load('ecv_colormap.mat');
ecv_color_map = a.c;
colormap(ecv_color_map)