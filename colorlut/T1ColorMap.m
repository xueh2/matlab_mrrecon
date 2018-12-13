function t1_color_map = T1ColorMap(use_cmap)

if(nargin<1)
    use_cmap = 1;
end

load t1_colormap
t1_color_map = c;

if(use_cmap)
    colormap(t1_color_map)
end