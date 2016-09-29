function [feFOV, peFOV, sliceThickness] = findFOVFromConfig(header);
% function [feFOV, peFOV] = findFOVFromConfig(header);
%
% function to find the FOV for fe and pe direction
%
% input:
%     header the headers.config part of [headers,protocol_header]=read_dat_headers(datfilename)
% output:
%     feFOV, peFOV the FOV in mm for fe and pe

str = '"ReadFoV">  { <Precision> ';
ind = strfind(header, str);
len = numel(str);
feFOV = str2num(header(ind+len+3:ind+len+12));

str = '"PhaseFoV">  { <Precision> ';
ind = strfind(header, str);
len = numel(str);
peFOV = str2num(header(ind+len+3:ind+len+12));

str = '<ParamArray."SliceThickness">';
ind = strfind(header, str);
len = numel(str);
sliceThickness = str2num(header(ind(1)+len+138:ind(1)+len+138+16));
