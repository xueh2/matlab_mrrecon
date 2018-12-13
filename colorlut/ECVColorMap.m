function ecv_color_map = ECVColorMap(use_cmap)

if(nargin<1)
    use_cmap = 1;
end

a = load('ecv_colormap.mat');
ecv_color_map = a.c;

if(use_cmap)
    colormap(ecv_color_map)
end