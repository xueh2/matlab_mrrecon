function [headers, imagePositionPatient2D, imageOrientationPatient2D] = CreateHeader2DFrameOf3DVolume(volume3D, header)
% compute header for every 2D frame of a 3D volume
% volume3D: [COL LIN SLC], 3D volume
% header: the MrFtk format header of 3D volume
% headers: the cell array storing the MrFtk format header for every 2D frame

N = size(volume3D, 3)
headers = cell(N, 1);
imagePositionPatient2D = zeros(N, 3);
imageOrientationPatient2D = zeros(N, 6);

imageOrientationPosition = [header.orientationPatient(1, :) header.orientationPatient(2, :)];
for i=1:N    
    header2D = header;
    header2D.sizeZ = 1;
    [wx, wy, wz] = Image2WorldMrFtk(header, 0, 0, i-1);    
    header2D = Dicom2HeaderMrFtk(volume3D(:, :, i), [header.spacingX header.spacingY header.spacingZ], [wx wy wz], imageOrientationPosition);
    
    imagePositionPatient2D(i, :) = [wx wy wz];
    imageOrientationPatient2D(i, :) = imageOrientationPosition;
    
    headers{i} = header2D;
end
