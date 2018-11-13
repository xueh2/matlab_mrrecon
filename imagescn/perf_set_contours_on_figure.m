
function fh = perf_set_contours_on_figure(fh, s, Seg, debugDir, contourDir)

%     h = openfig('T:\ReconResults\ML_Perf\template_roi.fig')
%     ROI_info_table = getappdata(h, 'ROI_Info_Table');

    ah = load('T:\ReconResults\ML_Perf\template_roi.mat');
    ROI_info_table = ah.ROI_info_table;
    
%     Seg = load(fullfile(contourDir, 'Seg'));

    endo0 = Seg.Seg(1).endo-1;
    endo1 = Seg.Seg(2).endo-1;
    endo2 = Seg.Seg(3).endo-1;

    epi0 = Seg.Seg(1).epi-1;
    epi1 = Seg.Seg(2).epi-1;
    epi2 = Seg.Seg(3).epi-1;
    
    endo0(:,1:2) = endo0(:,2:-1:1);
    endo1(:,1:2) = endo1(:,2:-1:1);
    endo2(:,1:2) = endo2(:,2:-1:1);
    
    epi0(:,1:2) = epi0(:,2:-1:1);
    epi1(:,1:2) = epi1(:,2:-1:1);
    epi2(:,1:2) = epi2(:,2:-1:1);
    
    ROI_info_table2 = ROI_info_table;

    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 1, 1, endo0+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 2, 1, epi0+1);
    
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 3, 2, endo1+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 4, 2, epi1+1);
    
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 5, 3, endo2+1);
    ROI_info_table2 = set_roi_info_table(ROI_info_table2, 6, 3, epi2+1);

    setappdata(fh, 'ROI_Info_Table', ROI_info_table2);

    ROI_info_table = ROI_info_table2;
    save(fullfile(contourDir, 'roi_adjust'), 'ROI_info_table');
    
    h_axes=flipud(findobj(fh,'type','axes'));

%     setappdata(h_axes(1),'ImageData', s(:,:,:,1));
%     h_image=findobj(h_axes(1),'type','image');
%     set(h_image,'cdata', s(:,:,1,1));
% 
%     setappdata(h_axes(2),'ImageData', s(:,:,:,2));
%     h_image=findobj(h_axes(2),'type','image');
%     set(h_image,'cdata', s(:,:,1,2));
% 
%     setappdata(h_axes(3),'ImageData', s(:,:,:,3));
%     h_image=findobj(h_axes(3),'type','image');
%     set(h_image,'cdata', s(:,:,1,3));

end

function rtable = set_roi_info_table(rtable, r, c, endo_epi)
    if(size(endo_epi,1)>2)
        rtable(r,c).ROI_Exists = 1;
        rtable(r,c).ROI_x_original = endo_epi(1:2:end,2);
        rtable(r,c).ROI_y_original = endo_epi(1:2:end,1);
        rtable(r,c).ROI_x_coordinates = endo_epi(1:2:end,2)';
        rtable(r,c).ROI_y_coordinates = endo_epi(1:2:end,1)';
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
