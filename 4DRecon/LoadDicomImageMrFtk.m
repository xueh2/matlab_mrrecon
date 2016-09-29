function [data, header] = LoadDicomImageMrFtk(filename)
% load a dicom image as MrFtk format
% data : [row col slice]

info = dicominfo(filename);
data = dicomread(filename);
header = CreateFtkHeaderInfo(data, [info.PixelSpacing' info.SliceThickness]);

% set the coordiante fields
header.positionPatient = info.ImagePositionPatient;

% row vector, 1st dimension, left-right
header.orientationPatient(1,:) = info.ImageOrientationPatient(1:3);

% col vector, 2nd dimension, up-down
header.orientationPatient(2,:) = info.ImageOrientationPatient(4:6);

% norm vector, 3nd dimension cross(r, c)
header.orientationPatient(3,:) = cross(header.orientationPatient(1,:), header.orientationPatient(2,:));
