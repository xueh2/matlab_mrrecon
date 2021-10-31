
function h = t1t2_plot_contours_on_images(t1, t2, t1w, t2w, pd, pts, endoC, epiC, t1_t2_values)
% h = t1t2_plot_contours_on_images(t1, t2, t1w, t2w, pd, pts, endoC, epiC)

    S = size(t1);

    if(numel(S)==2)
        SLC = 1;
    else
        SLC = S(3);
    end
    
    scale_factor = 16;

    h = figure;
    if(t1_t2_values(1)<1000)
        imagescn(t1, [0 1300], [1 SLC], scale_factor); T1ColorMap;
    else
        imagescn(t1, [0 2000], [1 SLC], scale_factor); T1ColorMap;
    end
    plot_contour_pt(h, pts, endoC, epiC, SLC);

    h_t2 = figure;
    imagescn(t2, [0 120], [1 SLC], scale_factor); T1ColorMap
    plot_contour_pt(h_t2, pts, endoC, epiC, SLC);

end

function plot_contour_pt(h, pts, endoC, epiC, SLC)
    h_axes=flipud(findobj(h,'type','axes'));

    contour_offset = 0;

    for slc=1:SLC
        axes(h_axes(slc))

        endo = endoC(:,:,slc);
        hold on    
        plot(endo(:,2)+contour_offset, endo(:,1)+contour_offset, 'b', 'LineWidth', 2);
        hold off

        epi = epiC(:,:,slc);
        hold on    
        plot(epi(:,2)+contour_offset, epi(:,1)+contour_offset, 'r', 'LineWidth', 2);
        hold off

        rvi_x = pts(1, 2, slc)+1;
        rvi_y = pts(1, 1, slc)+1;
        hold on    
        plot(rvi_y+contour_offset, rvi_x+contour_offset, 'g+', 'MarkerSize', 30, 'LineWidth', 3);
        hold off
    end
end