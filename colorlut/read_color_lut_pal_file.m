function [header, rgb]=read_color_lut_pal_file(filename);
% function [header, rgb]=read_color_lut_pal_file(filename);
%
% function to read the Siemens color LUT pal file (.pal)
%

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  NIH NHLBI                          *
%     ***************************************

fid=fopen(filename);

% read header (1st 5 lines)
header=[];
for i=1:5
    tline = fgetl(fid);
    header=strvcat(header,tline);
end

i=1;
while 1
    tline = fgetl(fid);
    if ~ischar(tline); fclose(fid); return; end
% line consists:
% [1] 0 , 0 , 145

    % get rgb values
    x=isspace(tline);
    ind1=findstr(tline,']'); % bracket
    if isempty(ind1); fclose(fid); return; end
    ind2=findstr(tline,','); % commas
    rgb(i,1)=str2num(deblank(tline(ind1(1)+1:ind2(1)-1)));
    rgb(i,2)=str2num(deblank(tline(ind2(1)+1:ind2(2)-1)));
    rgb(i,3)=str2num(deblank(tline(ind2(2)+1:end)));    

    % increment
    i=i+1;
end
fclose(fid);

return
