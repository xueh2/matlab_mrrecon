function ps_color_map = PSColorMap

a = load('ps_colormap.mat');
ps_color_map = a.PS_colormap;
colormap(ps_color_map)