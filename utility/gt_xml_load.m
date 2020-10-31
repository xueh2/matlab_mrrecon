function c = gt_xml_load(file)
% read in xml load
% c = gt_xml_load(file)

fid = fopen(file, 'r');

A = fread(fid);
A = char(A);
if(A(end)=='>')
else
    A = A(1:end-1);
end

fclose(fid);

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