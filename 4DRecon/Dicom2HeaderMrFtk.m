function header = Dicom2HeaderMrFtk(data, pixelSpacing, ImagePositionPatient, ImageOrientationPatient)
% Generate the MrFtk header in world coordinate as the dicom coordinate
% data : [row col slice]

header = CreateFtkHeaderInfo(data, pixelSpacing);

% set the coordiante fields
header.positionPatient = ImagePositionPatient;

% row vector, 1st dimension, left-right
header.orientationPatient(1,:) = ImageOrientationPatient(1:3);

% col vector, 2nd dimension, up-down
header.orientationPatient(2,:) = ImageOrientationPatient(4:6);

% norm vector, 3nd dimension cross(r, c)
header.orientationPatient(3,:) = cross(header.orientationPatient(1,:), header.orientationPatient(2,:));
