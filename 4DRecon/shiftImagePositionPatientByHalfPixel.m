
function imagePositionPatientShifted = shiftImagePositionPatientByHalfPixel(imagePositionPatient, imageOrientationPatient)
% here the patient position points to the upper-left corner of the first pixel, not the center of the first pixel,
% while MrFtk coordinates requires that image origin being the center of the first pixel
% thus, some correction is needed

imagePositionPatientShifted = imagePositionPatient;

N = size(imagePositionPatient, 1);

for i=1:N
    
    pos = imagePositionPatient(i, :);
    
    rowV = imageOrientationPatient(i, 1:3);
    colV = imageOrientationPatient(i, 4:6);
    
    posShifted = pos + 0.5*rowV + 0.5*colV;
    
    imagePositionPatientShifted(i, :) = posShifted;
end
