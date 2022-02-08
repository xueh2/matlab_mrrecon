function c = gt_xml_load_from_string(A)
% read in xml load from string
% c = gt_xml_load_from_string(file)

p = randi(100000);
suffix = datestr(now,'HHMMSSFFF');

if(isunix())
    fname = ['/tmp/temp_' num2str(p) '_' suffix '.attrib'];
else
    fname = ['c:/temp/temp_' num2str(p) '_' suffix  'attrib'];
end

fid = fopen(fname, 'w');
fprintf(fid, '%s', A);
fclose(fid);

c = xml2struct(fname);

delete(fname);