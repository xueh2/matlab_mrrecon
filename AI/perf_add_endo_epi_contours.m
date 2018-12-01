
function perf_add_endo_epi_contours(Gd, contourDir)

%     h = openfig('T:\ReconResults\ML_Perf\template_roi.fig')
%     ROI_info_table = getappdata(h, 'ROI_Info_Table');
% 
    ah = load('T:\ReconResults\ML_Perf\template_roi2.mat');
    ROI_info_table = ah.ROI_info_table;
    
%     figure; imagescn(Gd, [0 1.5], [1 3], [12], 3);
    
    try
        endo0 = analyze75read(fullfile(contourDir, 'endo_0'));
        endo1 = analyze75read(fullfile(contourDir, 'endo_1'));
        endo2 = analyze75read(fullfile(contourDir, 'endo_2'));

        epi0 = analyze75read(fullfile(contourDir, 'epi_0'));
        epi1 = analyze75read(fullfile(contourDir, 'epi_1'));
        epi2 = analyze75read(fullfile(contourDir, 'epi_2'));
        
        rv0 = analyze75read(fullfile(contourDir, 'endo_epi_rv_rvi_rv_0'));
        rv1 = analyze75read(fullfile(contourDir, 'endo_epi_rv_rvi_rv_1'));
        rv2 = analyze75read(fullfile(contourDir, 'endo_epi_rv_rvi_rv_2'));

        rvi0 = analyze75read(fullfile(contourDir, 'prob_endo_epi_rv_rvi_0'));
        rvi1 = analyze75read(fullfile(contourDir, 'prob_endo_epi_rv_rvi_1'));
        rvi2 = analyze75read(fullfile(contourDir, 'prob_endo_epi_rv_rvi_2'));

        rvi0 = rvi0(:,:,end);
        rvi1 = rvi1(:,:,end);
        rvi2 = rvi2(:,:,end);
        
        [v, ind] = max(rvi0(:)); [rvi0_s, rvi0_e] = CCMS_Contour(rvi0, v*0.9, 4, 0);
        [v, ind] = max(rvi0(:)); [rvi1_s, rvi1_e] = CCMS_Contour(rvi1, v*0.9, 4, 0);
        [v, ind] = max(rvi0(:)); [rvi2_s, rvi2_e] = CCMS_Contour(rvi2, v*0.9, 4, 0);
        
        ptN = size(rvi0_s,1);        
        rvi0 = zeros(2*ptN, 2);
        for pt=1:ptN
            rvi0( 2*(pt-1)+1, :) = rvi0_s(pt, 2:-1:1);
            rvi0( 2*(pt-1)+2, :) = rvi0_e(pt, 2:-1:1);
        end
        
        ptN = size(rvi1_s,1);        
        rvi1 = zeros(2*ptN, 2);
        for pt=1:ptN
            rvi1( 2*(pt-1)+1, :) = rvi1_s(pt, 2:-1:1);
            rvi1( 2*(pt-1)+2, :) = rvi1_e(pt, 2:-1:1);
        end
        
        ptN = size(rvi2_s,1);        
        rvi2 = zeros(2*ptN, 2);
        for pt=1:ptN
            rvi2( 2*(pt-1)+1, :) = rvi2_s(pt, 2:-1:1);
            rvi2( 2*(pt-1)+2, :) = rvi2_e(pt, 2:-1:1);
        end
        
    catch
        disp(['cannot load all contours : ' contourDir]);
        return
    end
    
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
            pt_list = [1:5:num_pt num_pt];
        else
            pt_list = [1:num_pt];
        end
        
        rtable(r,c).ROI_x_original = endo_epi(pt_list,2);
        rtable(r,c).ROI_y_original = endo_epi(pt_list,1);
        
        rtable(r,c).ROI_x_coordinates = endo_epi(pt_list,2)';
        rtable(r,c).ROI_y_coordinates = endo_epi(pt_list,1)';
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
