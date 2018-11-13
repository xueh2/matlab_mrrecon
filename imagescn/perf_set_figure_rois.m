
function h = perf_set_figure_rois(roi_file, data, template_file)

    h = openfig(template_file)
    ROI_info_table_fig = getappdata(h, 'ROI_Info_Table');

    rh = load('D:\data\template\roi.mat');
    ROI_info_table = rh.ROI_info_table;
       
    ah = load(roi_file);
    try
        ROI_info_table_loaded = ah.ROI_info_table;
    catch
        ROI_info_table_loaded = ah.roi.ROI_info_table;
    end
    
    endo0 = [ROI_info_table_loaded(1,1).ROI_x_original'; ROI_info_table_loaded(1,1).ROI_y_original']';
    endo1 = [ROI_info_table_loaded(3,2).ROI_x_original'; ROI_info_table_loaded(3,2).ROI_y_original']';
    endo2 = [ROI_info_table_loaded(5,3).ROI_x_original'; ROI_info_table_loaded(5,3).ROI_y_original']';
    
    epi0 = [ROI_info_table_loaded(2,1).ROI_x_original'; ROI_info_table_loaded(2,1).ROI_y_original']';
    epi1 = [ROI_info_table_loaded(4,2).ROI_x_original'; ROI_info_table_loaded(4,2).ROI_y_original']';
    epi2 = [ROI_info_table_loaded(6,3).ROI_x_original'; ROI_info_table_loaded(6,3).ROI_y_original']';
    
    ROI_info_table = set_roi_info_table(ROI_info_table, 1, 1, endo0);
    ROI_info_table = set_roi_info_table(ROI_info_table, 2, 1, epi0);
    
    ROI_info_table = set_roi_info_table(ROI_info_table, 3, 2, endo1);
    ROI_info_table = set_roi_info_table(ROI_info_table, 4, 2, epi1);
    
    ROI_info_table = set_roi_info_table(ROI_info_table, 5, 3, endo2);
    ROI_info_table = set_roi_info_table(ROI_info_table, 6, 3, epi2);

    sr = size(ROI_info_table);
    for c = 1:sr(2)
        for r = 1:sr(1)
            ROI_info_table(r, c).ROI_Data = ROI_info_table_fig(r, c).ROI_Data;
        end
    end
    
    setappdata(h, 'ROI_Info_Table', ROI_info_table);

    % update image
    h_axes=flipud(findobj(h,'type','axes'));

    if(size(data,4)>1)
        setappdata(h_axes(1),'ImageData', data(:,:,:,1));
        h_image=findobj(h_axes(1),'type','image');
        set(h_image,'cdata', data(:,:,1,1));

        setappdata(h_axes(2),'ImageData', data(:,:,:,2));
        h_image=findobj(h_axes(2),'type','image');
        set(h_image,'cdata', data(:,:,1,2));

        setappdata(h_axes(3),'ImageData', data(:,:,:,3));
        h_image=findobj(h_axes(3),'type','image');
        set(h_image,'cdata', data(:,:,1,3));
    else
        h_image=findobj(h_axes(1),'type','image');
        set(h_image,'cdata', data(:,:,1));

        h_image=findobj(h_axes(2),'type','image');
        set(h_image,'cdata', data(:,:,2));

        h_image=findobj(h_axes(3),'type','image');
        set(h_image,'cdata', data(:,:,3));
    end
    
    endo0 = [ROI_info_table(1,1).ROI_x_original'; ROI_info_table(1,1).ROI_y_original']';
    endo1 = [ROI_info_table(3,2).ROI_x_original'; ROI_info_table(3,2).ROI_y_original']';
    endo2 = [ROI_info_table(5,3).ROI_x_original'; ROI_info_table(5,3).ROI_y_original']';
    
    epi0 = [ROI_info_table(2,1).ROI_x_original'; ROI_info_table(2,1).ROI_y_original']';
    epi1 = [ROI_info_table(4,2).ROI_x_original'; ROI_info_table(4,2).ROI_y_original']';
    epi2 = [ROI_info_table(6,3).ROI_x_original'; ROI_info_table(6,3).ROI_y_original']';
    
    % update contours
    endo0_h = ROI_info_table(1,1).ROI_Data(1,1);
    epi0_h = ROI_info_table(2,1).ROI_Data(1,1);    
    set(endo0_h, 'XData', endo0(:,1), 'YData', endo0(:,2));
    set(epi0_h, 'XData', epi0(:,1), 'YData', epi0(:,2));
    
    endo1_h = ROI_info_table(3,2).ROI_Data(1,1);
    epi1_h = ROI_info_table(4,2).ROI_Data(1,1);    
    set(endo1_h, 'XData', endo1(:,1), 'YData', endo1(:,2));
    set(epi1_h, 'XData', epi1(:,1), 'YData', epi1(:,2));
    
    endo2_h = ROI_info_table(5,3).ROI_Data(1,1);
    epi2_h = ROI_info_table(6,3).ROI_Data(1,1);    
    set(endo2_h, 'XData', endo2(:,1), 'YData', endo2(:,2));
    set(epi2_h, 'XData', epi2(:,1), 'YData', epi2(:,2));
    drawnow;
end

function rtable = set_roi_info_table(rtable, r, c, endo_epi)
    if(size(endo_epi,1)>2)
        rtable(r,c).ROI_Exists = 1;
        rtable(r,c).ROI_x_original = endo_epi(:,1);
        rtable(r,c).ROI_y_original = endo_epi(:,2);
        rtable(r,c).ROI_x_coordinates = endo_epi(:,1)';
        rtable(r,c).ROI_y_coordinates = endo_epi(:,2)';
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
