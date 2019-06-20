function v = getAttribXMLField2(xmlContent, vname)
% v = getAttribXMLField(xmlContent, vname)

for ii=1:numel(xmlContent)        
    if(strcmp(xmlContent{ii}.name.Text, vname)==1)            
        for j=1:numel(xmlContent{ii}.value)
            v(j) = str2double(xmlContent{ii}.value{j}.Text);
        end            
    end        
end