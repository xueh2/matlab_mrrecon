
function [S, h] = cine_generate_seg_from_manual_roi(Cine, roi, phs, plotFlag, im_roi, scale_factor)
% contours are starting with index 1
% allow an extra upsampling ratio

if(nargin<4)
    plotFlag = 1;
end

if(nargin<5)
    im_roi = [];
    ro_s = -1;
    ro_e = -1;
    e1_s = -1;
    e1_e = -1;
end

if(nargin<6)
    scale_factor = 1;
end

S = struct('im', [], 'roi', [], ... 
        ...
        'endo', [], 'epi', [], 'rv', [], 'rvi', [], ...
        'endo_mask', [], 'epi_mask', [], 'rv_mask', [], 'myo_mask', [], 'rvi_mask', [], ...
        'endo_epi_mask', [], 'endo_epi_rvi_mask', [], 'endo_epi_rv_mask', [], 'endo_epi_rv_rvi_mask', []);
    
if(~isempty(im_roi))
    
    ro_s = im_roi(1);
    ro_e = im_roi(2);
    e1_s = im_roi(3);
    e1_e = im_roi(4);
    
    num_rois = size(roi, 1);
    num_slc = size(roi, 2);
    
    for slc=1:num_slc
        for r=1:num_rois
            
            roiC = roi(r, slc);
            if(roiC.ROI_Exists == 1)
                roiC.ROI_x_original = roiC.ROI_x_original + e1_s - 1;
                roiC.ROI_y_original = roiC.ROI_y_original + ro_s - 1;
                
                roiC.ROI_x_coordinates = roiC.ROI_x_coordinates + e1_s - 1;
                roiC.ROI_y_coordinates = roiC.ROI_y_coordinates + ro_s - 1;
                
                roi(r, slc) = roiC;
            end            
        end
    end    
end
    
S.im = squeeze(Cine(:,:,phs,:));
S.roi = roi;

RO = size(Cine, 1);
E1 = size(Cine, 2);
PHS = size(Cine, 3);
SLC = size(Cine, 4);

endo_mask = zeros(RO, E1, SLC);
epi_mask = endo_mask;
rvi_mask = endo_mask;

h = -1;

for slc = 1:SLC
    
    im = S.im(:,:,slc);
    
    %roiC = roi(2*(slc-1)+1, slc);
    endo_epi_C = find_rois(roi, slc);
    
    endoC = endo_epi_C{1};
    epiC = endo_epi_C{2};
    rviC = endo_epi_C{3};
    
    if(~isempty(endoC))
        pt_x = endoC.ROI_x_original;
        pt_y = endoC.ROI_y_original;
        
        pt_x = scale_factor * (pt_x-1) + 1;
        pt_y = scale_factor * (pt_y-1) + 1;
        
        endo_mask(:,:,slc) = single(roipoly(im, pt_x, pt_y));
    end

    if(~isempty(epiC))
        pt_x = epiC.ROI_x_original;
        pt_y = epiC.ROI_y_original;
        
        pt_x = scale_factor * (pt_x-1) + 1;
        pt_y = scale_factor * (pt_y-1) + 1;
        
        epi_mask(:,:,slc) = single(roipoly(im, pt_x, pt_y));
    end
    
    if(~isempty(rviC))
        pt_x = rviC.ROI_x_coordinates;
        pt_y = rviC.ROI_y_coordinates;
        
        pt_x = mean(pt_x);
        pt_y = mean(pt_y);
        
        pt_x = scale_factor * (pt_x-1) + 1;
        pt_y = scale_factor * (pt_y-1) + 1;
        
        %rvi_mask(:,:,slc) = single(roipoly(im, pt_x, pt_y));
        [X, Y] = meshgrid(round(pt_x)-2:round(pt_x)+2, round(pt_y)-2:round(pt_y)+2);
        rvi_mask(Y, X, slc) = 1;
    end
    
%     if(plotFlag)
%         h(slc) = plot_mask(im, endo_mask(:,:,slc), epi_mask(:,:,slc), [], rvi_mask(:,:,slc));
%     end
end

S.endo_mask = endo_mask;
S.epi_mask = epi_mask;
S.rvi_mask = rvi_mask;

if(plotFlag)
    h = plot_mask_all(S.im, S.endo_mask, S.epi_mask, S.rv_mask, S.rvi_mask, 12)
end

end

function h = plot_mask_all(cine, endo_mask, epi_mask, rv_mask, rvi_mask, scale_factor)

    SLC = size(cine, 3);
    
    h = figure;
    imagescn(cine, [1.1*min(cine(:)) 0.7*max(cine(:))], [3 ceil(SLC/3)], scale_factor);
    h_axes=flipud(findobj(h,'type','axes'));

    for slc=1:SLC    
        if(isempty(rv_mask))
            plot_mask(h_axes(slc), endo_mask(:,:,slc), epi_mask(:,:,slc), [], rvi_mask(:,:,slc));
        else
            plot_mask(h_axes(slc), endo_mask(:,:,slc), epi_mask(:,:,slc), rv_mask(:,:,slc), rvi_mask(:,:,slc));
        end
    end    
end

function plot_mask(ha, endo_mask, epi_mask, rv_mask, rvi_mask)
    axes(ha);
    [endo_s, endo_e] = CCMS_Contour(endo_mask, 0.5, 4, 0);
    [epi_s, epi_e] = CCMS_Contour(epi_mask, 0.5, 4, 0);

    if(~isempty(rv_mask))
        [rv_s, rv_e] = CCMS_Contour(rv_mask, 0.5, 4, 0);
    end
    if(~isempty(rvi_mask))
        [rvi_s, rvi_e] = CCMS_Contour(rvi_mask, 0.5, 4, 0);
    end
    
%     h = figure;imagescn(fmap);
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

function endo_epi_C = find_rois(roi, slc)

    endoC = [];
    epiC = [];
    rviC = [];

    if(isempty(roi))
        endo_epi_C = {endoC, epiC, rviC};
        return
    end
    
    num_rois = size(roi, 1);
       
    roi_ind = [];
    for ii=1:num_rois
        roiC = roi(ii, slc);
        if(roiC.ROI_Exists == 1)
            roi_ind = [roi_ind; ii];
        end
    end
    
    if(numel(roi_ind)==3)
        
        num_pts = [roi(roi_ind(1), slc).ROI_pixels, roi(roi_ind(2), slc).ROI_pixels, roi(roi_ind(3), slc).ROI_pixels];
        [pts, sorted_ind] = sort(num_pts);
        
        rviC = roi(roi_ind(sorted_ind(1)), slc);
        endoC = roi(roi_ind(sorted_ind(2)), slc);
        epiC = roi(roi_ind(sorted_ind(3)), slc);
        
%         if(roi(roi_ind(1), slc).ROI_pixels>roi(roi_ind(2), slc).ROI_pixels)
%             endoC = roi(roi_ind(2), slc);
%             epiC = roi(roi_ind(1), slc);
%         else
%             endoC = roi(roi_ind(1), slc);
%             epiC = roi(roi_ind(2), slc);
%         end   
%         
%         rviC = roi(roi_ind(3), slc);
    end
    
    if(numel(roi_ind)==2)
        
        if(roi(roi_ind(1), slc).ROI_pixels>roi(roi_ind(2), slc).ROI_pixels)
            endoC = roi(roi_ind(2), slc);
            epiC = roi(roi_ind(1), slc);
        else
            endoC = roi(roi_ind(1), slc);
            epiC = roi(roi_ind(2), slc);
        end        
    end
    
    if(numel(roi_ind)==1)        
        endoC = [];
        epiC = roi(roi_ind(1), slc);
    end
    
    endo_epi_C = {endoC, epiC, rviC};
end