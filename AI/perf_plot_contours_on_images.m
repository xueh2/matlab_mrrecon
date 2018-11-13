
function h = perf_plot_contours_on_images(Gd, endoC, epiC, rvC, rvi, windowing)

if (nargin<6)
    windowing = [0 1.5];
end

S = size(Gd);

scale_factor = 16;

if(S(3)>8)
    SLC = S(4);
    h = figure;
    imagescn(Gd, windowing, [1 SLC], scale_factor, 3);
else
    SLC = S(3);
    h = figure;
    imagescn(Gd, windowing, [1 SLC], scale_factor); PerfColorMap;
end

h_axes=flipud(findobj(h,'type','axes'));

for slc=1:SLC
    axes(h_axes(slc))
    
    if(numel(endoC)>=slc)
        endo = endoC{slc};
        hold on    
        plot(endo(:,2)+1, endo(:,1)+1, 'b', 'LineWidth', 2);
        hold off
    end
    
    if(numel(epiC)>=slc)
        epi = epiC{slc};
        hold on    
        plot(epi(:,2)+1, epi(:,1)+1, 'r', 'LineWidth', 2);
        hold off
    end
    
    if(numel(rvC)>=slc)
        rv = rvC{slc};
        hold on    
        plot(rv(:,2)+1, rv(:,1)+1, 'w', 'LineWidth', 2);
        hold off
    end 
    
    if(size(rvi, 1)>=slc)
        rvi_x = rvi(slc, 1);
        rvi_y = rvi(slc, 2);
        hold on    
        plot(rvi_y+1, rvi_x+1, 'm+', 'MarkerSize', 20);
        hold off
    end 
end