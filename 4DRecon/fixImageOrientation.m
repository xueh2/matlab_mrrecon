function [imagePositionPatientFixed, imageOrientationPatientFixed, dataFixed] = fixImageOrientation(imagePositionPatient, imageOrientationPatient, data, pixelSpacing)
% fix the image orientation if needed

imagePositionPatientFixed = imagePositionPatient;
imageOrientationPatientFixed = imageOrientationPatient;

numOfFrames = size(data, 3);

dataFixed = data;

cost = size(numOfFrames, 1);
cost(1) = 2.0;

frameCorrected = [];
for k=2:numOfFrames
    k
    posPre = imagePositionPatientFixed(k-1, :);
    rowVectorPre = imageOrientationPatientFixed(k-1, 1:3);
    colVectorPre = imageOrientationPatientFixed(k-1, 4:6);
    
    pos = imagePositionPatient(k, :);
    rowVector = imageOrientationPatient(k, 1:3);
    colVector = imageOrientationPatient(k, 4:6);
    
    cost(k) = dot(rowVector, rowVectorPre) + dot(colVectorPre, colVector);
    
    if ( cost(k) < 0.1 )
        % need correction
        frameCorrected = [frameCorrected k];        
        dataFrame = data(:,:,k);        
        [dataCorr, posCorr, rowCorr, colCorr] = performFrameOrientationCorrection(dataFrame, pixelSpacing, pos, rowVector,colVector, posPre, rowVectorPre,colVectorPre);
        
        dataFixed(:,:,k) = dataCorr;
        imagePositionPatientFixed(k, :) = posCorr;
        imageOrientationPatientFixed(k, :) = [rowCorr(1:3) colCorr(1:3)];
    end    
end