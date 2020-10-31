function attrib = gt_parse_image_attrib(xmlContent)
% parse the gt image attrib after reading in attrib file using xmlContent = gt_xml_load(attrib_file)
% attrib = gt_parse_image_attrib(xmlContent)

C = xmlContent.ismrmrdMeta.meta;

attrib = struct();

for n=1:numel(C)
    a_name = C{n}.name.Text;
    a_v = [];
    if(numel(C{n}.value)==1)
        if(~isempty(str2num(C{n}.value.Text)))
            a_v = [a_v str2num(C{n}.value.Text)];
        else
            a_v = [a_v {C{n}.value.Text}];
        end
    else
        for m=1:numel(C{n}.value)
            if(~isempty(str2num(C{n}.value{1}.Text)))
                a_v = [a_v str2num(C{n}.value{m}.Text)];
            else
                a_v = [a_v {C{n}.value{m}.Text}];
            end
        end
    end
    
    attrib = setfield(attrib, a_name, a_v);
end