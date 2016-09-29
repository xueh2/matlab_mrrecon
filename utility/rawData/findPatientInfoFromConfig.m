function [patientName, birthday, gender, patientID, studyID, seriesID, seqDescrip] = findPatientInfoFromConfig(header);
% function [patientName, birthday, gender] = findPatientInfoFromConfig(header);
%
% function to find the patientName, birthday, gender
%
% input:
%     header the headers.config part of [headers,protocol_header]=read_dat_headers(datfilename)
% output:
%     patientName, birthday, gender

str = '<ParamString."tPatientName">  { ';
ind = strfind(header, str);
len = numel(str);
patientName = header(ind+len+1:ind+len+16);

str = '<ParamString."PatientBirthDay">  { "';
ind = strfind(header, str);
len = numel(str);
birthday = header(ind(1)+len:ind(1)+len+7);

str = '<ParamLong."PatientSex">  { ';
ind = strfind(header, str);
len = numel(str);
gender = str2num(header(ind(1)+len:ind(1)+len+2));

str = '<ParamString."PatientLOID">  { "';
ind = strfind(header, str);
len = numel(str);
patientID = header(ind+len:ind+len+11);

str = '<ParamString."StudyLOID">  { "';
ind = strfind(header, str);
len = numel(str);
studyID = header(ind+len:ind+len+11);

str = '<ParamString."SeriesLOID">  { "';
ind = strfind(header, str);
len = numel(str);
seriesID = header(ind+len:ind+len+11);

str = '<ParamString."SequenceDescription">  { "';
ind = strfind(header, str);
len = numel(str);
seqDescrip = header(ind+len:ind+len+25);
