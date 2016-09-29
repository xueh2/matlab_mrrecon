

function data4D = LoadAndCorrectDicomImages(dataDir, plotFlag)
%
% this function load the data folder and sorted the images by acquisition time
% if needed, the image orientation will be corrected and dicom patient position and image orientation will be recomputed
%
% Inputs:
%    dataDir : folder storing the 4D datasets in dicom format (*.ima or *.dcm files)
%    plotFlag : if 1, then plot the information related to the datasets
%
% Output:
%    data4D: it is a structure including following fields:
%    all fields are stored after the sorting    
%
%    filenames : dicom file name
%    info : dicom info loaded
%    sliceLocation
%    triggerTime
%    acquisitionTime 
%    dataAcquired : loaded image data before correction
%    imagePositionPatient : position patient before correction
%    imageOrientationPatient : orientation patient before correction
%    imagePositionPatientFixed  : position patient after correction
%    imageOrientationPatientFixed  : orientation patient after correction
%    dataAcquiredFixed : image data after correction
%
%     ***************************************
%     *  Hui Xue (hui-xue@siemens.com       *
%     *  2012-04                            *
%     ***************************************

%% get the folder information
cd(dataDir)
currDir = dataDir;

[names, numFile] = findFILE(currDir, '*.ima');
if ( numFile == 0 )
    [names, numFile] = findFILE(currDir, '*.dcm');
end

if ( numFile == 0 )
    error(['No images are found at ' currDir]);
end

disp('-------------------------------------------------------');

%% load the dicom images
disp('Loading the images ... ');

info = dicominfo(names{1});
data = dicomread(names{1});
header = CreateFtkHeaderInfo(data, [info.PixelSpacing' info.SliceThickness]);

imagePositionPatient = zeros(numFile, 3);
imageOrientationPatient = zeros(numFile, 6);
sliceLocation = zeros(numFile, 1); % in the unit of mm
triggerTime = zeros(numFile, 1); % in the unit of ms
acquisitionTime = zeros(numFile, 1);  % in the unit of ms
dataAcquired = zeros(header.sizeY, header.sizeX, numFile);
infoes = cell(numFile, 1);
imageNumber = zeros(numFile, 1);
filenames = cell(numFile, 1);

tic

% load image info
for k=1:numFile
    [pathstr, name, ext] = fileparts(names{k});
    infoes{k} = dicominfo(names{k});
    imageNumber(k) = infoes{k}.InstanceNumber;
    timeInSeconds = ConvertDicomAcquisitionTime2Seconds(info.AcquisitionTime);
    acquisitionTime(k) = timeInSeconds;
end

% sort by imagenumber
[acquisitionTime2, ind] = sort(acquisitionTime);

imageNumber2 = imageNumber(ind);

% load image content
for k=1:numFile    
    [pathstr, name, ext] = fileparts(names{ind(k)});
    info = infoes{ind(k)};
    sliceLocation(k) = info.SliceLocation;
    triggerTime(k) = info.TriggerTime;
    timeInSeconds = ConvertDicomAcquisitionTime2Seconds(info.AcquisitionTime);
    acquisitionTime(k) = timeInSeconds;
    imagePositionPatient(k, :) = info.ImagePositionPatient;
    imageOrientationPatient(k, :) = info.ImageOrientationPatient;
    filenames{k} = names{ind(k)};
    
    disp([num2str(ind(k)) ' - ' name ' - ' num2str(sliceLocation(k)) ' - ' num2str(triggerTime(k)) ' - ' info.AcquisitionTime ' - ' num2str(timeInSeconds)]);
    
    dataX = dicomread(names{ind(k)});    
    dataAcquired(:,:,k) = dataX;
end

headerAll = CreateFtkHeaderInfo(dataAcquired, [info.PixelSpacing' info.SliceThickness]);
Matlab_SaveAnalyze(single(dataAcquired), headerAll, fullfile(currDir, 'dataAcquired.hdr'));
disp(['Loading images : ' num2str(toc)]);

%% here the patient position points to the upper-left corner of the first pixel, not the center of the first pixel,
% while MrFtk coordinates requires that image origin being the center of the first pixel
% thus, some correction is needed
% imagePositionPatientShifted = shiftImagePositionPatientByHalfPixel(imagePositionPatient, imageOrientationPatient);

%% check and fix the data properties
disp('Correcting the image orientation ... ');
pixelSpacing = [info.PixelSpacing' info.SliceThickness];
[imagePositionPatientFixed, imageOrientationPatientFixed, dataAcquiredFixed] = fixImageOrientation(imagePositionPatient, imageOrientationPatient, dataAcquired, pixelSpacing);

save dataReadyForProcessing filenames info sliceLocation triggerTime acquisitionTime imagePositionPatient ... 
    imageOrientationPatient imagePositionPatientFixed imageOrientationPatientFixed dataAcquiredFixed dataAcquired headerAll

headerAll = CreateFtkHeaderInfo(dataAcquiredFixed, [info.PixelSpacing' info.SliceThickness]);
Matlab_SaveAnalyze(single(dataAcquiredFixed), headerAll, fullfile(currDir, 'dataAcquiredFixed.hdr'));

% if necessary, correct acquisition time
acquisitionTime = acquisitionTime - acquisitionTime(1);
acquisitionDelta = acquisitionTime(2:end) - acquisitionTime(1:end-1);
ind = find(acquisitionDelta>1)
if ( ~isempty(ind) )
    
    numOfGaps = numel(ind);
    acquisitionTime2 = acquisitionTime;
    
    for g=1:numOfGaps
        gap = acquisitionTime(ind(g)+1) - acquisitionTime(ind(g))
        gap2 = acquisitionTime(ind(g)) - acquisitionTime(ind(g)-1)
        acquisitionTime2(ind(g)+1:end) = acquisitionTime2(ind(g)+1:end) - gap + gap2;
    end
    
    figure;
    scatter(1:numFile, acquisitionTime2, '+');
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Image number');
    ylabel('Acquisition time corrected');
    
    acquisitionTimeOri = acquisitionTime;
    acquisitionTime = acquisitionTime2;
end

data4D = struct('filenames', filenames, 'info', info, 'sliceLocation', sliceLocation, ... 
    'triggerTime', triggerTime, 'acquisitionTime', acquisitionTime, 'dataAcquired', dataAcquired, ... 
    'imagePositionPatient', imagePositionPatient, 'imageOrientationPatient', imageOrientationPatient, ...
    'imagePositionPatientFixed', imagePositionPatientFixed, 'imageOrientationPatientFixed', imageOrientationPatientFixed, ...
    'dataAcquiredFixed', dataAcquiredFixed);

if ( plotFlag )
    figure;
    scatter(triggerTime, sliceLocation, '.');
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Trigger Time - ms');
    ylabel('Slice Location - mm');

    figure;
    scatter(imageNumber2, triggerTime, '.');
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Image number');
    ylabel('Trigger time - ms');

    figure;
    scatter(imageNumber2, sliceLocation, '.');
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Image number');
    ylabel('Slice location - mm');

    figure;
    scatter(imageNumber2, acquisitionTime, '.');
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Image number');
    ylabel('Acquisition time');

    figure;
    scatter(sliceLocation, acquisitionTime, '.');
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Slice Location - mm');
    ylabel('Acquisition Time - ms');
    
    figure;
    hold on
    scatter3(imagePositionPatientFixed(:, 1), imagePositionPatientFixed(:,2), imagePositionPatientFixed(:,3));
    hold off
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('ImagePositionPatientX - mm');
    ylabel('ImagePositionPatientY - mm');
    zlabel('ImagePositionPatientZ - mm');
    
    figure;
    hold on
    scatter3(imageOrientationPatient(:, 1), imageOrientationPatient(:,2), imageOrientationPatient(:,3));
    hold off
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('imageOrientationPatient row vector x - mm');
    ylabel('imageOrientationPatient row vector y - mm');
    zlabel('imageOrientationPatient row vector z - mm');
    
    plotStep = floor(numFile/100);
    if ( plotStep < 5 ) plotStep = 5; end
    ha = -1;
    aviobj = avifile(fullfile(currDir, 'slice3D.avi'));
    for k=1:plotStep:numFile
        [d, h] = LoadDicomImageMrFtk(filenames{k});
        ha = plotMrFtkImageIn3D(d, h, ha, 200, 400);
        frame = getframe(gcf);
        aviobj = addframe(aviobj,frame);
    end
    aviobj = close(aviobj);
    
    plotStep = floor(numFile/100);
    if ( plotStep < 5 ) plotStep = 5; end
    ha = -1;
    aviobj2 = avifile(fullfile(currDir, 'slice3DFixedOrientation.avi'));
    for k=1:plotStep:numFile
        h = Dicom2HeaderMrFtk(data, pixelSpacing, imagePositionPatientFixed(k,:), imageOrientationPatientFixed(k,:));
        d = dataAcquiredFixed(:,:,k);
        ha = plotMrFtkImageIn3D(d, h, ha, 200, 400);
        frame = getframe(gcf);
        aviobj2 = addframe(aviobj2,frame);
    end
    aviobj2 = close(aviobj2);    
end
