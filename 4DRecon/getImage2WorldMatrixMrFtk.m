function m = getImage2WorldMatrixMrFtk(header)
% image to world matrix for MrFtk

translateM = eye(4,4);
translateM(1,4) = header.positionPatient(1);
translateM(2,4) = header.positionPatient(2);
translateM(3,4) = header.positionPatient(3);

rotateM = eye(4,4);
rotateM(1:3,1:3) = header.orientationPatient;
rotateM = rotateM';

scaleM = eye(4,4);
scaleM(1,1) = header.spacingX;
scaleM(2,2) = header.spacingY;
scaleM(3,3) = header.spacingZ;

m = translateM*rotateM*scaleM;

