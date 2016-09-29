function [age, weight] = findPatientInfoFromDicom(header);
% function [age] = findPatientInfoFromDicom(header);
%
% function to find the age
%
% input:
%     header the headers.Dicom part of [headers,protocol_header]=read_dat_headers(datfilename)
% output:
%     patientName, birthday, gender

str = '<ParamDouble."flPatientAge">  { <Precision> 6  ';
ind = strfind(header, str);
len = numel(str);
age = str2num(header(ind+len:ind+len+9));

str = '<ParamDouble."flUsedPatientWeight">  { <Precision> 6  ';
ind = strfind(header, str);
len = numel(str);
weight = str2num(header(ind+len:ind+len+9));
