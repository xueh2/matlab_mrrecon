function [dataSlc, headerSlc] = getSLC(data4D, header, s)
% data4D: [COL LIN SLC PHS]
% s starts from 1

headerSlc = header;
dataSlc = squeeze(data4D(:,:,s,:));
headerSlc.sizeZ = size(dataSlc, 3);
headerSlc.sizeT = 1;

% compute the new image position patient
[wx, wy, wz ] = Image2WorldMrFtk(header, 0, 0, s-1);
headerSlc.positionPatient = [wx wy wz];   
