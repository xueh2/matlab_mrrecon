
function cine_add_endo_epi_contours(Cine, endo, epi, phs, contourDir)

%     h = openfig('T:\ReconResults\ML_Perf\template_roi.fig')
%     ROI_info_table = getappdata(h, 'ROI_Info_Table');
% 

    if(isunix())
        baseDir = '/mnt/Lab-Kellman/'
    else
        baseDir = '\\hl-share\RawMRI\Lab-Kellman\'
    end

    RO = size(Cine,1);
    E1 = size(Cine,2);
    PHS = size(Cine,3);
    SLC = size(Cine,4);

    ah = load(fullfile(baseDir, 'RawData', 'MachinLearning_Labelled_data', 'Cine', 'roi.mat'));
    ROI_info_table = ah.ROI_info_table;
        
%     figure; imagescn(Cine, [], [], [12], 3);
    
    ROI_info_table2 = ROI_info_table;

    if (SLC>size(ROI_info_table, 2))
        SLC = size(ROI_info_table,2);
    end
    
    for slc = 1:SLC
        if(~isempty(endo{slc, phs}))
            ROI_info_table2 = set_roi_info_table(ROI_info_table2, 2*(slc-1)+1, slc, endo{slc, phs}+1);
        else
            ROI_info_table2 = set_roi_info_table(ROI_info_table2, 2*(slc-1)+1, slc, []);
        end
        
        if(~isempty(epi{slc, phs}))
            ROI_info_table2 = set_roi_info_table(ROI_info_table2, 2*(slc-1)+2, slc, epi{slc, phs}+1);
        else
            ROI_info_table2 = set_roi_info_table(ROI_info_table2, 2*(slc-1)+2, slc, []);
        end
    end    
    
    % setappdata(h, 'ROI_Info_Table', ROI_info_table2);

    ROI_info_table = ROI_info_table2(1:2*SLC, :);
    save(fullfile(contourDir, ['roi_phs' num2str(phs)]), 'ROI_info_table');
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
