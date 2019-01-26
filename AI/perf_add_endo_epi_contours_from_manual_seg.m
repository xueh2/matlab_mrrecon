
function perf_add_endo_epi_contours_from_manual_seg(contourDir, endo, epi, rv, rvi)
% perf_add_endo_epi_contours_from_manual_seg(contourDir, endo, epi, rv, rvi)
% input endo/epi/rv contours are 1 based cell array
% rvi: [SLC 2] rv insertion point

    ah = load('T:\ReconResults\ML_Perf\template_roi2.mat');
    ROI_info_table = ah.ROI_info_table;
    
%     figure; imagescn(Gd, [0 1.5], [1 3], [12], 3);
    
    endo0 = endo{1}-1;
    endo1 = endo{2}-1;
    endo2 = endo{3}-1;

    epi0 = epi{1}-1;
    epi1 = epi{2}-1;
    epi2 = epi{3}-1;

    if(~isempty(rv{1}))
        rv0 = rv{1}-1;
    else
        rv0 = endo0+5;
    end
    
    if(~isempty(rv{2}))
        rv1 = rv{2}-1;
    else
        rv1 = endo1+5;
    end
    
    if(~isempty(rv{3}))
        rv2 = rv{3}-1;
    else
        rv2 = endo2+5;
    end
    
    rvi0 = rvi(1,:)-1;
    rvi1 = rvi(2,:)-1;
    rvi2 = rvi(3,:)-1;
    
    ROI_info_table2 = ROI_info_table;

    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 1, 1, endo0+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 2, 1, epi0+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 3, 1, rv0+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 4, 1, rvi0+1);
    
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 5, 2, endo1+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 6, 2, epi1+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 7, 2, rv1+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 8, 2, rvi1+1);
    
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 9,  3, endo2+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 10, 3, epi2+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 11, 3, rv2+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 12, 3, rvi2+1);

    % setappdata(h, 'ROI_Info_Table', ROI_info_table2);

    ROI_info_table = ROI_info_table2;
    save(fullfile(contourDir, 'roi'), 'ROI_info_table');
%     h_axes=flipud(findobj(h,'type','axes'));
% 
%     s1 = Gd(:,:,:,1);
%     s2 = Gd(:,:,:,2);
%     s3 = Gd(:,:,:,3);
%     
%     setappdata(h_axes(1),'ImageData', s1);
%     h_image=findobj(h_axes(1),'type','image');
%     set(h_image,'cdata', s1(:,:,1));
% 
%     setappdata(h_axes(2),'ImageData', s2);
%     h_image=findobj(h_axes(2),'type','image');
%     set(h_image,'cdata', s2(:,:,1));
% 
%     setappdata(h_axes(3),'ImageData', s3);
%     h_image=findobj(h_axes(3),'type','image');
%     set(h_image,'cdata', s3(:,:,1));

end

function rtable = set_roi_info_table(rtable, r, c, endo_epi)
    if(size(endo_epi,1)>2)
        rtable(r,c).ROI_Exists = 1;
        
        num_pt = size(endo_epi, 1);
        if(num_pt>100)
            stepSize = round(num_pt/80);
            pt_list = [1:stepSize:num_pt num_pt];
        elseif(num_pt>50)
            pt_list = [1:5:num_pt num_pt];
        else
            pt_list = [1:num_pt];
        end
        
        rtable(r,c).ROI_x_original = endo_epi(pt_list,1);
        rtable(r,c).ROI_y_original = endo_epi(pt_list,2);
        
        rtable(r,c).ROI_x_coordinates = endo_epi(pt_list,1)';
        rtable(r,c).ROI_y_coordinates = endo_epi(pt_list,2)';
%     else
%         rtable(r,c).ROI_Exists = 0;
%         rtable(r,c).ROI_Info = [];
%         rtable(r,c).ROI_Data = [];
%         rtable(r,c).ROI_x_original = [];
%         rtable(r,c).ROI_y_original = [];
%         rtable(r,c).ROI_x_coordinates = [];
%         rtable(r,c).ROI_y_coordinates = [];        
    end
end
