function All_ROIs = plot_ai_cine_seg_barts_mosaic_with_rois(sax_mosaiced_array, prefix, ind)
% plot_ai_cine_seg_barts_mosaic_with_rois(sax_mosaiced_array, prefix)

if (nargin<3)
    ind = [];
end

PHS = size(sax_mosaiced_array, 3);

% load ROI
All_ROIs = [];
if(~isempty(ind))
    fname = [prefix num2str(ind) '.attrib'];
    xmlcontent = gt_xml_load(fname);
    All_ROIs = getXMLField(xmlcontent, 'GT_ROI');    
else    
    for n=1:PHS
        fname = [prefix num2str(n-1) '.attrib'];
        xmlcontent = gt_xml_load(fname);
        ROIs = getXMLField(xmlcontent, 'GT_ROI');
        
        All_ROIs = [All_ROIs ROIs];
    end
end

if(~isempty(ind))    
    im = sax_mosaiced_array(:,:,ind+1);
    h = figure; imagescn(im, [min(im(:)) 0.25*max(im(:))], [], 12);
    hold on
    for n=1:numel(All_ROIs)
        C = All_ROIs{n};
        color = C(1:3);
        num_pt = (numel(C)-4)/2;
        ro = C(5:2:end); e1 = C(6:2:end);
        plot(e1+1, ro+1, 'LineWidth', C(4), 'Color', color);
    end
    hold off    
end

end

function ROIs = getXMLField(xmlContent, vname)

    ROIs = [];
    C = xmlContent.ismrmrdMeta.meta;
    for ii=1:numel(C)        
        if(~isempty(strfind(C{ii}.name.Text, vname)))
            v = [];
            for j=1:numel(C{ii}.value)
                v(j) = str2double(C{ii}.value{j}.Text);
            end    
            ROIs = [ROIs {v}];
        end        
    end
end