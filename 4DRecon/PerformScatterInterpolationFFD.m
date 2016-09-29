function [volumeBFFD, headerBFFD] = PerformScatterInterpolationFFD(binData, pixelSpacing, ...
    imagePositionPatient, imageOrientationPosition, ...
    dstPixelSpacing, numOfSubdivision, tmpFolder, backgroundThres)
% perform the scatter interpolation using FFD
% binData : [COL LIN SLC PHS], reconed bin data
% image position and orientation for every 2D frame
% pixelSpacing : the resoltuion for binData
% dstPixelSpacing : the required resoltuion for recon volume

% scatter interpolation

if ( nargin < 8 )
    backgroundThres = -1;
end

% find the recon size from the first phase
data = binData(:,:,:,1);
header = CreateFtkHeaderInfo(data, pixelSpacing);
N = size(data, 3)
volume = cell(N, 1);
headers = cell(N, 1);
for i=1:N    
    volume{i} = data(:,:,i);    
    header2D = header;
    header2D.sizeZ = 1;    
    header2D = Dicom2HeaderMrFtk(volume{i}, pixelSpacing, imagePositionPatient(i, :), imageOrientationPosition(i, :));           
    headers{i} = header2D;
end
[boundingBox, headerResampled] = findBoundingBoxes(volume, headers, dstPixelSpacing);

% ensure the maximal level of approximation
controlPointSpacing = 0.5*[headerResampled.spacingX*header.sizeX  headerResampled.spacingY*header.sizeY  headerResampled.spacingZ*header.sizeZ];

% recon every phase
NPhs = size(binData, 4)
volumeBFFD = zeros(headerResampled.sizeY, headerResampled.sizeX, headerResampled.sizeZ, NPhs);

for phs=1:NPhs
    disp(['Processing phase ' num2str(phs) ' ... ']);
    data = binData(:,:,:,phs);
    
    header = CreateFtkHeaderInfo(data, pixelSpacing);

    N = size(data, 3)

    volume = cell(N, 1);
    headers = cell(N, 1);
    
    usedN = 1;
    for i=1:N    
        
        value = double(data(:,:,i));
        
        minV = min(value(:));
        maxV = max(value(:));
        
        if ( maxV - minV < 1 )
            continue;
        end
        
        volume{usedN} = value;
        header2D = header;
        header2D.sizeZ = 1;    
        header2D = Dicom2HeaderMrFtk(volume{usedN}, pixelSpacing, imagePositionPatient(i, :), imageOrientationPosition(i, :));           
        headers{usedN} = header2D;
        
        usedN = usedN + 1;
    end
    
    if ( N ~= usedN-1 )
        disp('Black image found ... ');
    end
    
    volume = volume(1:usedN-1);
    headers = headers(1:usedN-1);
    
    tic; volumeResampled = Matlab_ScatterInterpolationFFD(volume, headers, headerResampled, controlPointSpacing, numOfSubdivision, backgroundThres, tmpFolder); toc
      
    if ( ~isempty(tmpFolder) )
        fileName = fullfile(tmpFolder, ['volumeResampled' '_' num2str(phs) '.hdr']);
        Matlab_SaveAnalyze(single(volumeResampled), headerResampled, fileName);
    end
    
    volumeBFFD(:,:,:,phs) = volumeResampled;
end

fileName = fullfile(tmpFolder, ['volumeResampled' '4D' '.hdr']);
headerBFFD = headerResampled;
headerBFFD.sizeT = size(volumeBFFD, 4);
Matlab_SaveAnalyze(single(volumeBFFD), headerBFFD, fileName);
