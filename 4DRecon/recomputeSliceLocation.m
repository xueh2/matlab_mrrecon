function sliceLocationNew = recomputeSliceLocation(data, pixelSpacing, sliceLocation, imagePositionPatient, imageOrientationPatient)
% Recompute the slice location if the orignal slice location is not monotolic
% if the orignal slice location is not monotolic, then the rotating sampling strategy is used

sliceLocationNew = sliceLocation;

minSlc = min(sliceLocation(:));
maxSlc = max(sliceLocation(:));

isRotating = 0;
if ( (sliceLocation(1)~=minSlc) & (sliceLocation(1)~=maxSlc) )
    isRotating = 1;
end

if ( (sliceLocation(end)~=minSlc) & (sliceLocation(end)~=maxSlc) )
    isRotating = 1;
end

if ( isRotating )
    
    N = numel(sliceLocation);
    sliceLocationNew(1) = 0;
    
    header = Dicom2HeaderMrFtk(data, pixelSpacing, imagePositionPatient(1,:), imageOrientationPatient(1,:));
    [wc1, wc2, wc3, wc4] = computeFourCorner(data, header);
    
    sliceVectorFirst = wc4-wc1;
    sliceVectorFirst = sliceVectorFirst / norm(sliceVectorFirst);
    for k=2:N       
        header = Dicom2HeaderMrFtk(data, pixelSpacing, imagePositionPatient(k,:), imageOrientationPatient(k,:));
        [wc1, wc2, wc3, wc4] = computeFourCorner(data, header);
        sliceVector = wc4-wc1;
        sliceVector = sliceVector / norm(sliceVector);        
        v = dot(sliceVectorFirst, sliceVector);
        sliceLocationNew(k) = acos(v)*180/pi;
    end    
end

