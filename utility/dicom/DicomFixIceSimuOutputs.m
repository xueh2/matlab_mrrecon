function startSeriesNum = DicomFixIceSimuOutputs(dicomDir, ExampleDicomImg, startSeriesNum)
% fix the ICE outputs of dicom images to enable the import to the scanner

% gather all info
info = dicominfo(ExampleDicomImg);

info.PatientName
info.PatientID
info.PatientBirthDate
info.PatientAge
info.PatientSex
info.PatientSize
info.PatientWeight

info.StudyDate
info.StudyDescription
info.StudyID
info.StudyInstanceUID

info.ReferencedImageSequence
info.ReferringPhysicianName
info.RequestedProcedureDescription

info.ImplementationVersionName

info.AccessionNumber

[names, num] = findFILE(dicomDir, '*.dcm');

if ( num == 0 )
    return;
end

seriesNum = [];

for ii=1:num
    disp(names{ii});
    
    info2 = dicominfo(names{ii});
    data = dicomread(info2);
    
    % fix the dicom header
    info2.PatientName = info.PatientName;
    info2.PatientID = info.PatientID;
    info2.PatientBirthDate = info.PatientBirthDate;
    info2.PatientAge = info.PatientAge;
    info2.PatientSex = info.PatientSex;
    info2.PatientSize = info.PatientSize;
    info2.PatientWeight = info.PatientWeight;

    info2.StudyDate = info.StudyDate;
    info2.StudyDescription = info.StudyDescription;
    info2.StudyID = info.StudyID;
    info2.StudyInstanceUID = info.StudyInstanceUID;

    info2.ReferencedImageSequence = info.ReferencedImageSequence;
    info2.ReferringPhysicianName = info.ReferringPhysicianName;
    info2.RequestedProcedureDescription = info.RequestedProcedureDescription;

    info2.ImplementationVersionName = info.ImplementationVersionName;

    info2.AccessionNumber = info.AccessionNumber;

    info2.SeriesNumber = startSeriesNum + info2.SeriesNumber;
    
    if ( isfield(info2, 'SeriesDescription') )
        info2.SeriesDescription = [info2.SeriesDescription '_ICESIMU'];
    else
        info2.SeriesDescription = [info2.ProtocolName '_GT_ICESIMU'];
    end
    
    if ( isfield(info2, 'ImageComments') )
        info2.ImageComments = [info2.ImageComments '_ICESIMU'];
    else
        info2.ImageComments = [info2.ProtocolName '_GT_ICESIMU'];
    end
    
    seriesNum = [seriesNum info2.SeriesNumber];
    
    % save the dicom images
    dicomwrite(data, names{ii}, info2);
end

startSeriesNum = max(seriesNum(:)) + 1;

