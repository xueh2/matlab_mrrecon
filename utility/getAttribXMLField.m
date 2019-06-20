function v = getAttribXMLField(xmlContent, vname)
% v = getAttribXMLField(xmlContent, vname)

for ii=1:numel(xmlContent)        
    if(strcmp(xmlContent(ii).meta(1).name, vname)==1)            
        for j=1:numel(xmlContent(ii).meta)
            v(j) = str2double(xmlContent(ii).meta(j).value);
        end            
    end        
end