
%% test the scatter interpolation using FFD

cd D:\cc_views\hxue_pccs000608ws_snapshot_view2\SCR_MR_CV\MrFtk\data\RegData

[data, header] = Matlab_LoadAnalyze('3D_target_rigid_F.hdr');
header.orientationPatient = eye(3);

pixelSpacing = [header.spacingX header.spacingY header.spacingZ];
% load a fack dicom info
info = dicominfo('D:\data\4DRecon_20111206\20120223\selectedSeries\series32\dicom\2CV_VOLUNTEER2.MR.PLM-AW_CARD_PESPE.0032.3590.2012.02.23.19.14.59.375000.15880902.IMA');
header2 = Dicom2HeaderMrFtk(data, pixelSpacing, info.ImagePositionPatient, info.ImageOrientationPatient);

header.positionPatient = header2.positionPatient;
header.orientationPatient = header2.orientationPatient;

N = header.sizeZ;

volume = cell(N, 1);
headers = cell(N, 1);

for i=1:N
    
    volume{i} = data(:,:,i);
    
    header2D = header;
    header2D.sizeZ = 1;
    
    [wx, wy, wz ] = Image2WorldMrFtk(header, 0, 0, i-1);
    header2D.positionPatient = [wx wy wz];
    
    headers{i} = header2D;
end

plotStep = 1;
ha = -1;
for k=1:plotStep:N
    d = volume{k};
    h = headers{k}
    ha = plotMrFtkImageIn3D(d, h, ha, 1024, 1024);
end
        
headerDst = header;
controlPointSpacing = 16*[header.spacingX  header.spacingY  header.spacingZ];
numOfSubdivision = 6;
tmpFolder = 'd:/temp/'

cd D:\cc_views\hxue_pccs000608ws_snapshot_view2\SCR_MR_CV\MrFtk\prod\bin\vc10_Debug\Debug

volumeResampled = Matlab_ScatterInterpolationFFD(volume, headers, headerDst, controlPointSpacing, numOfSubdivision, tmpFolder);

Matlab_SaveAnalyze(single(volumeResampled), headerDst, fullfile(tmpFolder, 'volumeResampled.hdr'));

%% do the real test
cd D:\data\4DRecon_20111206\20120223\selectedSeries\series32

load HeartBeatPickingResults

% plotStep = 1;
% saveAvi = 1;
% aviName = 'd:/temp/3DRadialVolume.avi';
% 
% pixelSpacing = [1.822916626930200   1.822916626930200   7.000000000000000]
% 
% slc = 25;
% Render4DVolume(binAveAll, pixelSpacing, bestHeartBeatsImagePositionPatient, bestHeartBeatsImageOrientationPatient, slc, -1, plotStep, 200, 400, 1/20);
% 
% phs = 6;
% Render4DVolume(binAveAll, pixelSpacing, bestHeartBeatsImagePositionPatient, bestHeartBeatsImageOrientationPatient, -1, phs, plotStep, 200, 400, 1/20, saveAvi, aviName);
% 
% % recon phs 6
% data = binAveAll(:,:,:,6);
% size(data)
% 
% header = CreateFtkHeaderInfo(data, pixelSpacing);
% 
% N = size(data, 3)
% 
% volume = cell(N, 1);
% headers = cell(N, 1);
% for i=1:N    
%     volume{i} = data(:,:,i);    
%     header2D = header;
%     header2D.sizeZ = 1;    
%     header2D = Dicom2HeaderMrFtk(volume{i}, pixelSpacing, bestHeartBeatsImagePositionPatient(i, :), bestHeartBeatsImageOrientationPatient(i, :));           
%     headers{i} = header2D;
% end
% 
% % find the bounding box of all slices
% [boundingBox, headerResampled] = findBoundingBoxes(volume, headers, [1.822916626930200 1.822916626930200 1.822916626930200]);
% 
% % draw the slices
plotStep = 1;
ha = -1;
for k=1:plotStep:N
    d = volume{k};
    h = headers{k}
    ha = plotMrFtkImageIn3D(d, h, ha, 1024, 1024);
end
% 
% volumeResampled = zeros(headerResampled.sizeY, headerResampled.sizeX, headerResampled.sizeZ);
% 
% % draw the bounding box
% h = plotBoundingBoxOf3DVolume(volumeResampled, headerResampled, ha);

% scatter interpolation

controlPointSpacing = 0.5*[headerResampled.spacingX*header.sizeX  headerResampled.spacingY*header.sizeY  headerResampled.spacingZ*header.sizeZ];
numOfSubdivision = 9;
tmpFolder = 'd:/temp/'

NPhs = size(binAveAll, 4)

volumeBFFD = zeros(headerResampled.sizeY, headerResampled.sizeX, headerResampled.sizeZ, NPhs)

for phs=1:NPhs
    disp(['Processing phase ' num2str(phs) ' ... ']);
    data = binAveAll(:,:,:,phs);
    
    header = CreateFtkHeaderInfo(data, pixelSpacing);

    N = size(data, 3)

    volume = cell(N, 1);
    headers = cell(N, 1);
    for i=1:N    
        volume{i} = data(:,:,i);    
        header2D = header;
        header2D.sizeZ = 1;    
        header2D = Dicom2HeaderMrFtk(volume{i}, pixelSpacing, bestHeartBeatsImagePositionPatient(i, :), bestHeartBeatsImageOrientationPatient(i, :));           
        headers{i} = header2D;
    end
    
    [boundingBox, headerResampled] = findBoundingBoxes(volume, headers, [1.822916626930200 1.822916626930200 1.822916626930200]);
    volumeResampled = zeros(headerResampled.sizeY, headerResampled.sizeX, headerResampled.sizeZ);

    volumeResampled = Matlab_ScatterInterpolationFFD(volume, headers, headerResampled, controlPointSpacing, numOfSubdivision, tmpFolder);
      
    fileName = fullfile(tmpFolder, ['volumeResampled' '_' num2str(phs) '.hdr']);
    Matlab_SaveAnalyze(single(volumeResampled), headerResampled, fileName);
    
    volumeBFFD(:,:,:,phs) = volumeResampled;
end


for phs=1:NPhs
    disp(['Processing phase ' num2str(phs) ' ... ']);
      
    fileName = fullfile(tmpFolder, ['volumeResampled' '_' num2str(phs) '.hdr']);
    [volumeResampled, header] = Matlab_LoadAnalyze(fileName);
    
    volumeBFFD(:,:,:,phs) = volumeResampled;
end

fileName = fullfile(tmpFolder, ['volumeResampled' '4D' '.hdr']);
header4D = headerResampled;
header4D.sizeT = size(volumeBFFD, 4);
Matlab_SaveAnalyze(single(volumeBFFD), header4D, fileName);

save4DVolume(volumeBFFD, headerResampled, tmpFolder, 'FFD', 1, 200, 400, 1/24);
