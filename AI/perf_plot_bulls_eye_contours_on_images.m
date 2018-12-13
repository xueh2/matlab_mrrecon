
function [h, h2] = perf_plot_bulls_eye_contours_on_images(Gd, bulls_eye, windowing)

    if (nargin<3)
        windowing = [0 1.5];
    end

    S = size(Gd);

    scale_factor = 16;

    if(S(3)>8)
        SLC = S(4);
        h = figure;
        imagescn(Gd, windowing, [1 SLC], scale_factor, 3);
        
        h2 = figure;
        imagescn(Gd, windowing, [1 SLC], scale_factor, 3);
    else
        SLC = S(3);
        h = figure;
        imagescn(Gd, windowing, [1 SLC], scale_factor); PerfColorMap;
        
        h2 = figure;
        imagescn(Gd, windowing, [1 SLC], scale_factor); PerfColorMap;
    end

    contour_offset = 0;
    
    h_axes=flipud(findobj(h,'type','axes'));
    for slc=1:SLC

        slc_plot = bulls_eye(slc).slice_index;

        axes(h_axes(slc_plot))

        plot_c(bulls_eye(slc).endo_contours, 'r', contour_offset);
        plot_c(bulls_eye(slc).epi_contours, 'g', contour_offset);
    end

    h_axes=flipud(findobj(h2,'type','axes'));
    for slc=1:SLC

        slc_plot = bulls_eye(slc).slice_index;

        axes(h_axes(slc_plot))

        plot_c(bulls_eye(slc).contours, 'b', contour_offset);
    end
end

function plot_c(endo_all, c, contour_offset)
    hold on    
    for ii=1:numel(endo_all)
        endo = endo_all{ii};
        plot(endo(:,1)+contour_offset, endo(:,2)+contour_offset, c, 'LineWidth', 2); 
        text(mean(endo(:,1)), mean(endo(:,2)), [num2str(ii) ')'], 'FontSize', 16, 'Color', [1 1 1]);
    end
    hold off
end