function posVec = ComputeDicomPosVecFromICE(posVecICE, rowVec, colVec, pixelSpacingRow, pixelSpacingCol, noOfRows, noOfCols)
% compute the dicom position vector from ICE
% posVec = ComputeDicomPosVecFromICE(posVecICE, rowVec, ColVec, pixelSpacingRow, pixelSpacingCol, noOfRows, noOfCols)

aHalfRowFoV = pixelSpacingCol * noOfCols / 2;
aHalfColFoV = pixelSpacingRow * noOfRows / 2;

posVec = posVecICE - rowVec * aHalfRowFoV - colVec * aHalfColFoV;

