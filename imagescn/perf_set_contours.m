
function h = perf_set_contours_on_figure(debugDir, contourDir)

    h = openfig('T:\ReconResults\ML_Perf\template_roi.fig')
    ROI_info_table = getappdata(h, 'ROI_Info_Table');

    ah = load('T:\ReconResults\ML_Perf\template_roi.mat');
    ROI_info_table = ah.ROI_info_table;
    
    s1 = analyze75read(fullfile(debugDir, 'perf_moco_upsampled_for_seg_0'));
    s2 = analyze75read(fullfile(debugDir, 'perf_moco_upsampled_for_seg_1'));
    s3 = analyze75read(fullfile(debugDir, 'perf_moco_upsampled_for_seg_2'));

    figure; imagescn(cat(4, s1, s2, s3), [0 1.5], [1 3], [12], 3);
    
    fmap0 = analyze75read(fullfile(debugDir, 'flow_maps_for_seg_0'));
    fmap1 = analyze75read(fullfile(debugDir, 'flow_maps_for_seg_1'));
    fmap2 = analyze75read(fullfile(debugDir, 'flow_maps_for_seg_2'));

    try
        endo0 = analyze75read(fullfile(contourDir, 'endo_0'));
        endo1 = analyze75read(fullfile(contourDir, 'endo_1'));
        endo2 = analyze75read(fullfile(contourDir, 'endo_2'));

        epi0 = analyze75read(fullfile(contourDir, 'epi_0'));
        epi1 = analyze75read(fullfile(contourDir, 'epi_1'));
        epi2 = analyze75read(fullfile(contourDir, 'epi_2'));
    catch
        Seg = load(fullfile(contourDir, 'Seg'));
        
        endo0 = Seg.Seg(1).endo-1;
        endo1 = Seg.Seg(2).endo-1;
        endo2 = Seg.Seg(3).endo-1;
        
        epi0 = Seg.Seg(1).epi-1;
        epi1 = Seg.Seg(2).epi-1;
        epi2 = Seg.Seg(3).epi-1;
    end
    
    ROI_info_table2 = ROI_info_table;

    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 1, 1, endo0+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 2, 1, epi0+1);
    
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 3, 2, endo1+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 4, 2, epi1+1);
    
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 5, 3, endo2+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 6, 3, epi2+1);

    setappdata(h, 'ROI_Info_Table', ROI_info_table2);

    ROI_info_table = ROI_info_table2;
    save(fullfile(contourDir, 'roi'), 'ROI_info_table');
    
    h_axes=flipud(findobj(h,'type','axes'));

    setappdata(h_axes(1),'ImageData', s1);
    h_image=findobj(h_axes(1),'type','image');
    set(h_image,'cdata', s1(:,:,1));

    setappdata(h_axes(2),'ImageData', s2);
    h_image=findobj(h_axes(2),'type','image');
    set(h_image,'cdata', s2(:,:,1));

    setappdata(h_axes(3),'ImageData', s3);
    h_image=findobj(h_axes(3),'type','image');
    set(h_image,'cdata', s3(:,:,1));

end

function rtable = set_roi_info_table(rtable, r, c, endo_epi)
    if(size(endo_epi,1)>2)
        rtable(r,c).ROI_Exists = 1;
        
        num_pt = size(endo_epi, 1);
        pt_list = [1:5:num_pt num_pt];
        
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
