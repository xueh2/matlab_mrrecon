function [volumeResampled, headerResampled] = resampleVolumeForFigure(raw, header)

d = zeros(header.sizeY, header.sizeX, header.sizeZ+2);
d(:,:,2:header.sizeZ+1) = raw;
d(:,:,1) = raw(:,:,1);
d(:,:,end) = raw(:,:,end);

header2 = header;
header2.sizeZ = size(d, 3);
header2.positionPatient = [0 0 -8];

dstVolume = zeros([4*header2.sizeY 4*header2.sizeX header.sizeZ+4]);
headerDst = header2;
headerDst.sizeX = size(dstVolume, 2);
headerDst.sizeY = size(dstVolume, 1);
headerDst.sizeZ = size(dstVolume, 3);;
headerDst.spacingX = headerDst.spacingX/4;
headerDst.spacingY = headerDst.spacingY/4;

headerDst.positionPatient = [0 0 -16];

[volumeResampled, headerResampled] = resampleVolume(double(d), header2, dstVolume, headerDst);
volumeResampled = volumeResampled(:,:, 3:header.sizeZ+2);
            