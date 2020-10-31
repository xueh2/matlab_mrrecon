
function S = generate_seg_from_manual_roi(data, roi, plotFlag)
% contours are starting with index 1
% data: [RO E1 SLC/PHS]
% roi: contours for endo/epi/rvi ...
% S = generate_seg_from_manual_roi(data, roi, plotFlag)

S = struct('im', [], 'endo', [], 'epi', [], 'rv', [], 'rvi', [], ...
        'endo_mask', [], 'epi_mask', [], 'rv_mask', [], 'myo_mask', [], 'rvi_mask', [], ...
        'endo_epi_mask', [], 'endo_epi_rvi_mask', [], 'endo_epi_rv_mask', [], 'endo_epi_rv_rvi_mask', []);
    
S.im = data;

RO = size(data, 1);
E1 = size(data, 2);
SLC = size(data, 3);

S.endo_mask = zeros(RO, E1, SLC) -1;
S.epi_mask = zeros(RO, E1, SLC)-1;
S.rv_mask = zeros(RO, E1, SLC)-1;
S.myo_mask = zeros(RO, E1, SLC)-1;
S.rvi_mask = zeros(RO, E1, SLC)-1;
S.endo_epi_mask = zeros(RO, E1, SLC)-1;
S.endo_epi_rvi_mask = zeros(RO, E1, SLC)-1;
S.endo_epi_rv_mask = zeros(RO, E1, SLC)-1;
S.endo_epi_rv_rvi_mask = zeros(RO, E1, SLC)-1;
 
num_rois = size(roi, 1);

for slc=1:SLC
    
    valid_rois = [];
    valid_rois_pixels = [];
    for r=1:num_rois
        if(roi(r, slc).ROI_Exists)
            valid_rois = [valid_rois r];
            valid_rois_pixels = [valid_rois_pixels roi(r, slc).ROI_pixels];
        end
    end    
    
    endoC = [];
    epiC = [];
    rvC = [];
    rvi = [];
    
    if(numel(valid_rois)==2) % endo, epi
        if(valid_rois_pixels(1)>valid_rois_pixels(2))
            epiC = [roi(valid_rois(1), slc).ROI_x_original roi(valid_rois(1), slc).ROI_y_original];
            endoC = [roi(valid_rois(2), slc).ROI_x_original roi(valid_rois(2), slc).ROI_y_original];
        else
            epiC = [roi(valid_rois(2), slc).ROI_x_original roi(valid_rois(2), slc).ROI_y_original];
            endoC = [roi(valid_rois(1), slc).ROI_x_original roi(valid_rois(1), slc).ROI_y_original];
        end
    end
    
    if(numel(valid_rois)==3) % endo, epi, rvi
        [valid_rois_pixels_sorted, ind] = sort(valid_rois_pixels);
        epiC = [roi(valid_rois(ind(3)), slc).ROI_x_original roi(valid_rois(ind(3)), slc).ROI_y_original];
        endoC = [roi(valid_rois(ind(2)), slc).ROI_x_original roi(valid_rois(ind(2)), slc).ROI_y_original];
        rvi = [mean(roi(valid_rois(ind(1)), slc).ROI_x_original) mean(roi(valid_rois(ind(1)), slc).ROI_y_original)];
    end
    
    if(numel(valid_rois)>=2)
        [endo_mask, epi_mask, rv_mask, myo_mask, rvi_mask, ...
            endo_epi_mask, endo_epi_rv_mask, endo_epi_rvi_mask, endo_epi_rv_rvi_mask] = create_mask(endoC, epiC, rvC, rvi, data(:,:,slc));

        S.endo{slc} = endoC;
        S.epi{slc} = epiC;
        S.rv{slc} = rvC;
        S.rvi{slc} = rvi;

        S.endo_mask(:,:,slc) = endo_mask;
        S.epi_mask(:,:,slc) = epi_mask;

        if(~isempty(rv_mask))
            S.rv_mask(:,:,slc) = rv_mask;
        end

        S.rvi_mask(:,:,slc) = rvi_mask;

        S.myo_mask(:,:,slc) = myo_mask;
        S.endo_epi_mask(:,:,slc) = endo_epi_mask;
        S.endo_epi_rv_mask(:,:,slc) = endo_epi_rv_mask;
        S.endo_epi_rvi_mask(:,:,slc) = endo_epi_rvi_mask;
        S.endo_epi_rv_rvi_mask(:,:,slc) = endo_epi_rv_rvi_mask;
    end
end

if(plotFlag)
    for slc=1:SLC
        plot_mask(data(:,:,slc), S.endo_mask(:,:,slc), S.epi_mask(:,:,slc), S.rv_mask(:,:,slc), S.rvi_mask(:,:,slc));
    end
end

end

function [endo_mask, epi_mask, rv_mask, myo_mask, rvi_mask, endo_epi_mask, endo_epi_rv_mask, endo_epi_rvi_mask, endo_epi_rv_rvi_mask] = create_mask(endo, epi, rv, rvi, fmap)

    endo_mask = single(roipoly(fmap, endo(:,1), endo(:,2)));    
    epi_mask = single(roipoly(fmap, epi(:,1), epi(:,2)));   
    
    rv_mask = [];
    if(~isempty(rv))
        rv_mask = single(roipoly(fmap, rv(:,1), rv(:,2)));    
    end
    
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
            rv_i_x_all(x+(y-1)*3) = rv_i_x+(x-2);
            rv_i_y_all(x+(y-1)*3) = rv_i_y+(y-2);
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
    
    if(~isempty(rv_mask))
        ind = find(rv_mask(:)>0);
        endo_epi_rv_mask(ind) = 3;
    end
    
    endo_epi_rv_rvi_mask = endo_epi_rv_mask;
    endo_epi_rv_rvi_mask(rv_i_y_all, rv_i_x_all) = 4;
end

function plot_mask(data, endo_mask, epi_mask, rv_mask, rvi_mask)
    [endo_s, endo_e] = CCMS_Contour(endo_mask, 0.5, 4, 0);
    [epi_s, epi_e] = CCMS_Contour(epi_mask, 0.5, 4, 0);

    if(sum(rv_mask(:))>0)
        [rv_s, rv_e] = CCMS_Contour(rv_mask, 0.5, 4, 0);
    end
    if(sum(rvi_mask(:))>0)
        [rvi_s, rvi_e] = CCMS_Contour(rvi_mask, 0.5, 4, 0);
    end
    
    figure;imagescn(data, [], [], [4]);
    hold on
    for tt=1:size(endo_s,1)-1
        line([endo_s(tt,1) endo_e(tt,1)],[endo_s(tt,2) endo_e(tt,2)], ...
            'LineWidth', 1, 'Color', 'r'); 
    end
    for tt=1:size(epi_s,1)-1
        line([epi_s(tt,1) epi_e(tt,1)],[epi_s(tt,2) epi_e(tt,2)], ...
        'LineWidth', 1, 'Color', 'b'); 
    end 

    if(sum(rv_mask(:))>0)
        for tt=1:size(rv_s,1)-1
            line([rv_s(tt,1) rv_e(tt,1)],[rv_s(tt,2) rv_e(tt,2)], ...
            'LineWidth', 1, 'Color', 'm'); 
        end
    end
    if(sum(rvi_mask(:))>0)
        for tt=1:size(rvi_s,1)-1
            line([rvi_s(tt,1) rvi_e(tt,1)],[rvi_s(tt,2) rvi_e(tt,2)], ...
            'LineWidth', 3, 'Color', 'g'); 
        end
    end
    hold off
end