function [dataPhs, headerPhs] = getPHS(data4D, header, phs)
% data4D: [COL LIN SLC PHS]
% s starts from 1

headerPhs = header;
dataPhs = squeeze(data4D(:,:,:,phs));
headerPhs.sizeZ = size(dataPhs, 3);
headerPhs.sizeT = 1;
