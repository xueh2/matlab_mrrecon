function ps_color_map = PSColorMap(use_cmap)

if(nargin<1)
    use_cmap = 1;
end

a = load('PS_colormap.mat');
ps_color_map = a.PS_colormap;

if(use_cmap)
    colormap(ps_color_map)
end