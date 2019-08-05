
function ROI_info_table = gt_AI_add_contours_to_roi(data, C, roi_file)
% gt_AI_add_contours_to_roi(data, C)
% data: [RO E1], an image
% C: a cell array for every contour
% e.g. roi_file = 'T:\Share\MATLAB\MrRecon\AI\roi_2D.mat'

    RO = size(data,1);
    E1 = size(data,2);

    ah = load(roi_file);
    ROI_info_table = ah.ROI_info_table;
           
    ROI_info_table2 = ROI_info_table;

    num_c = numel(C);
    
    for r=1:num_c
    
        if(~isempty(C{r}))
            ROI_info_table2 = set_roi_info_table(ROI_info_table2, r, 1, C{r});
        else
            ROI_info_table2 = set_roi_info_table(ROI_info_table2, r, 1, []);
        end
    end
    
    ROI_info_table = ROI_info_table2(1:num_c, :);
end

function rtable = set_roi_info_table(rtable, r, c, endo_epi)
    if(~isempty(endo_epi) & size(endo_epi,1)>2)
        rtable(r,c).ROI_Exists = 1;
        
        num_pt = size(endo_epi, 1);
        if(num_pt>100)
            pt_list = [1:5:num_pt num_pt];
        else
            pt_list = [1:num_pt];
        end
        
        rtable(r,c).ROI_x_original = endo_epi(pt_list,1);
        rtable(r,c).ROI_y_original = endo_epi(pt_list,2);
        
        rtable(r,c).ROI_x_coordinates = endo_epi(pt_list,1)';
        rtable(r,c).ROI_y_coordinates = endo_epi(pt_list,2)';
    else
        rtable(r,c).ROI_Exists = 0;
        rtable(r,c).ROI_Info = [];
        rtable(r,c).ROI_Data = [];
        rtable(r,c).ROI_x_original = [];
        rtable(r,c).ROI_y_original = [];
        rtable(r,c).ROI_x_coordinates = [];
        rtable(r,c).ROI_y_coordinates = [];        
    end
end
