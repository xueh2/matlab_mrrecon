function ps_color_map = PSColorMap

a = load('PS_colormap.mat');
ps_color_map = a.PS_colormap;
colormap(ps_color_map)