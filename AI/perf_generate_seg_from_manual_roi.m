
function S = perf_generate_seg_from_manual_roi(roi, contour_roi, fmap, fmap_resized, fmap_resized_training, startROI, slc, plotFlag)
% contours are starting with index 1

S = struct('im', [], 'roi', [], ... 
        ...
        'endo_resized_training', [], 'epi_resized_training', [], 'rv_resized_training', [], 'rvi_resized_training', [], ...
        'endo_resized_training_mask', [], 'epi_resized_training_mask', [], 'rv_resized_training_mask', [], 'myo_resized_training_mask', [], 'rvi_resized_training_mask', [], ... 
        'endo_epi_resized_training_mask', [], 'endo_epi_rvi_resized_training_mask', [], 'endo_epi_rv_resized_training_mask', [], 'endo_epi_rv_rvi_resized_training_mask', [],  ...
        ...
        'endo_resized', [], 'epi_resized', [], 'rv_resized', [], 'rvi_resized', [], ...
        'endo_resized_mask', [], 'epi_resized_mask', [], 'rv_resized_mask', [], 'myo_resized_mask', [], 'rvi_resized_mask', [], ...
        'endo_epi_resized_mask', [], 'endo_epi_rvi_resized_mask', [], 'endo_epi_rv_resized_mask', [], 'endo_epi_rv_rvi_resized_mask', [],  ...
        ...
        'endo', [], 'epi', [], 'rv', [], 'rvi', [], ...
        'endo_mask', [], 'epi_mask', [], 'rv_mask', [], 'myo_mask', [], 'rvi_mask', [], ...
        'endo_epi_mask', [], 'endo_epi_rvi_mask', [], 'endo_epi_rv_mask', [], 'endo_epi_rv_rvi_mask', []);
    
S.im = fmap_resized_training;
S.roi = contour_roi;

ps_ro = contour_roi(1);
pe_ro = contour_roi(2);

ps_e1 = contour_roi(3);
pe_e1 = contour_roi(4);

ps_ro = ps_ro -1;
ps_e1 = ps_e1 -1;

rs = size(fmap_resized, 1)/size(fmap,1);

% --------------------------------------------

S.endo_resized_training = [roi.ROI_info_table(startROI,slc).ROI_x_original roi.ROI_info_table(startROI,slc).ROI_y_original];
S.epi_resized_training = [roi.ROI_info_table(startROI+1,slc).ROI_x_original roi.ROI_info_table(startROI+1,slc).ROI_y_original];
S.rv_resized_training = [roi.ROI_info_table(startROI+2,slc).ROI_x_original roi.ROI_info_table(startROI+2,slc).ROI_y_original];
S.rvi_resized_training = [ mean(roi.ROI_info_table(startROI+3,slc).ROI_x_original) mean(roi.ROI_info_table(startROI+3,slc).ROI_y_original)];

[S.endo_resized_training_mask, S.epi_resized_training_mask, S.rv_resized_training_mask, S.myo_resized_training_mask, S.rvi_resized_training_mask, S.endo_epi_resized_training_mask, S.endo_epi_rv_resized_training_mask, S.endo_epi_rvi_resized_training_mask, S.endo_epi_rv_rvi_resized_training_mask] = create_mask(S.endo_resized_training, S.epi_resized_training, S.rv_resized_training, S.rvi_resized_training, fmap_resized_training);

% --------------------------------------------

S.endo_resized = S.endo_resized_training;
S.endo_resized(:,1) = S.endo_resized_training(:,1) + ps_e1;
S.endo_resized(:,2) = S.endo_resized_training(:,2) + ps_ro;

S.epi_resized = S.epi_resized_training;
S.epi_resized(:,1) = S.epi_resized_training(:,1) + ps_e1;
S.epi_resized(:,2) = S.epi_resized_training(:,2) + ps_ro;

S.rv_resized = S.rv_resized_training;
S.rv_resized(:,1) = S.rv_resized_training(:,1) + ps_e1;
S.rv_resized(:,2) = S.rv_resized_training(:,2) + ps_ro;

S.rvi_resized = S.rvi_resized_training;
S.rvi_resized(1) = S.rvi_resized(1) + ps_e1;
S.rvi_resized(2) = S.rvi_resized(2) + ps_ro;

[S.endo_resized_mask, S.epi_resized_mask, S.rv_resized_mask, S.myo_resized_mask, S.rvi_resized_mask, S.endo_epi_resized_mask, S.endo_epi_rv_resized_mask, S.endo_epi_rvi_resized_mask, S.endo_epi_rv_rvi_resized_mask] = create_mask(S.endo_resized, S.epi_resized, S.rv_resized, S.rvi_resized, fmap_resized);

% --------------------------------------------

S.endo = remove_scaling(S.endo_resized, rs);
S.epi = remove_scaling(S.epi_resized, rs);
S.rv = remove_scaling(S.rv_resized, rs);
S.rvi = remove_scaling(S.rvi_resized, rs);
[S.endo_mask, S.epi_mask, S.rv_mask, S.myo_mask, S.rvi_mask, S.endo_epi_mask, S.endo_epi_rv_mask, S.endo_epi_rvi_mask, S.endo_epi_rv_rvi_mask] = create_mask(S.endo, S.epi, S.rv, S.rvi, fmap);

if(plotFlag)
    figure; imagescn(cat(3, S.endo_mask, S.epi_mask, S.rv_mask, S.myo_mask, S.rvi_mask, S.endo_epi_mask, S.endo_epi_rv_mask, S.endo_epi_rvi_mask, S.endo_epi_rv_rvi_mask), [], [3 3]);
    
    plot_mask(fmap, S.endo_mask, S.epi_mask, S.rv_mask, S.rvi_mask);
    plot_mask(fmap_resized, S.endo_resized_mask, S.epi_resized_mask, S.rv_resized_mask, S.rvi_resized_mask);
    plot_mask(fmap_resized_training, S.endo_resized_training_mask, S.epi_resized_training_mask, S.rv_resized_training_mask, S.rvi_resized_training_mask);
end

end

function [endo_mask, epi_mask, rv_mask, myo_mask, rvi_mask, endo_epi_mask, endo_epi_rv_mask, endo_epi_rvi_mask, endo_epi_rv_rvi_mask] = create_mask(endo, epi, rv, rvi, fmap)

    endo_mask = single(roipoly(fmap, endo(:,1), endo(:,2)));    
    epi_mask = single(roipoly(fmap, epi(:,1), epi(:,2)));    
    rv_mask = single(roipoly(fmap, rv(:,1), rv(:,2)));    

    myo_mask = single(epi_mask - endo_mask);
    ind = find(myo_mask(:)<0);
    myo_mask(ind) = 0;

    endo_epi_mask = endo_mask;
    endo_epi_mask(:) = 0;
    ind = find(endo_mask(:)>0);
    endo_epi_mask(ind) = 1;
    ind = find(myo_mask(:)>0);
    endo_epi_mask(ind) = 2;

    rv_i_x = round(rvi(1));
    rv_i_y = round(rvi(2));

    rv_i_x_all = zeros(9,1);
    rv_i_y_all = zeros(9,1);

    for y=1:3
        for x=1:3
            rv_i_x_all(x+(y-1)*3) = rv_i_x-x-2;
            rv_i_y_all(x+(y-1)*3) = rv_i_y-y-2;
        end
    end

    rvi_mask = endo_mask;
    rvi_mask(:) = 0;
    rvi_mask(rv_i_y_all, rv_i_x_all) = 1;

    endo_epi_rvi_mask = endo_epi_mask;
    endo_epi_rvi_mask(rv_i_y_all, rv_i_x_all) = 3;

    endo_epi_rv_mask = endo_mask;
    endo_epi_rv_mask(:) = 0;
    ind = find(endo_mask(:)>0);
    endo_epi_rv_mask(ind) = 1;
    ind = find(myo_mask(:)>0);
    endo_epi_rv_mask(ind) = 2;
    ind = find(rv_mask(:)>0);
    endo_epi_rv_mask(ind) = 3;

    endo_epi_rv_rvi_mask = endo_epi_rv_mask;
    endo_epi_rv_rvi_mask(rv_i_y_all, rv_i_x_all) = 4;
end

function endo = remove_scaling(endo_resized, rs)
    endo = (endo_resized-1)/rs + 1;
end

function plot_mask(fmap, endo_mask, epi_mask, rv_mask, rvi_mask)
    [endo_s, endo_e] = CCMS_Contour(endo_mask, 0.5, 4, 0);
    [epi_s, epi_e] = CCMS_Contour(epi_mask, 0.5, 4, 0);

    if(~isempty(rv_mask))
        [rv_s, rv_e] = CCMS_Contour(rv_mask, 0.5, 4, 0);
    end
    if(~isempty(rvi_mask))
        [rvi_s, rvi_e] = CCMS_Contour(rvi_mask, 0.5, 4, 0);
    end
    
    figure;imagescn(fmap, [0 8]); PerfColorMap;
    hold on
    for tt=1:size(endo_s,1)-1
        line([endo_s(tt,1) endo_e(tt,1)],[endo_s(tt,2) endo_e(tt,2)], ...
            'LineWidth', 1, 'Color', 'r'); 
    end
    for tt=1:size(epi_s,1)-1
        line([epi_s(tt,1) epi_e(tt,1)],[epi_s(tt,2) epi_e(tt,2)], ...
        'LineWidth', 1, 'Color', 'b'); 
    end 

    if(~isempty(rv_mask))
        for tt=1:size(rv_s,1)-1
            line([rv_s(tt,1) rv_e(tt,1)],[rv_s(tt,2) rv_e(tt,2)], ...
            'LineWidth', 1, 'Color', 'm'); 
        end
    end
    if(~isempty(rvi_mask))
        for tt=1:size(rvi_s,1)-1
            line([rvi_s(tt,1) rvi_e(tt,1)],[rvi_s(tt,2) rvi_e(tt,2)], ...
            'LineWidth', 3, 'Color', 'g'); 
        end
    end
    hold off
end