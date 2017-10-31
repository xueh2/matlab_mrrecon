function [headers, data] = DicomConverter2(home, dstDir)

disp('parsing the home directory ...')

subdirs = [];
subdirs = FindDeepSubDirs(home, subdirs, [])

% subdirs = [subdirs {home}];

info = struct('filename', cell(0), 'PatientName', '', 'SeriesNumber', -1, 'SliceLocation', -1,...
    'AcquisitionNumber', -1, 'ImageSize', [-1 -1], 'PixelSpacing', [0 0 0], 'InversionTime', -1);
fileInfo = [];

num = numel(subdirs);
for k=1:num
    
    D = dir( fullfile( subdirs{k} ,'*.*') );
    fileNum = numel(D);
    
    disp(['Read images at ' subdirs{k}])
    
    for s = 1:fileNum
        if ( D(s).isdir == 0 )
            
            currentFileName = fullfile(subdirs{k}, D(s).name)
            if ( isdicom(currentFileName) )
                info = struct('filename', '', 'PatientName', '', 'SeriesNumber', -1, 'SliceLocation', -1,...
                    'AcquisitionNumber', -1, 'ImageSize', [-1 -1], 'PixelSpacing', [0 0 0], 'InversionTime', -1);

                cinfo = dicominfo(currentFileName);
                info.filename = currentFileName;
                
                if ( isfield(cinfo, 'PatientName') )
                    info.PatientName = cinfo.PatientName.FamilyName;
                else
                    info.PatientName = 'Unknowun';
                end
                
                if ( isfield(cinfo, 'SeriesNumber') )
                    info.SeriesNumber = cinfo.SeriesNumber;
                end
                
                if ( isfield(cinfo, 'SliceLocation') )
                    info.SliceLocation = cinfo.SliceLocation;
                end
                
                if ( isfield(cinfo, 'AcquisitionNumber') )
                    info.AcquisitionNumber = cinfo.AcquisitionNumber;                
                end

                if ( isfield(cinfo, 'ImageNumber') )
                    info.AcquisitionNumber = cinfo.ImageNumber;                
                end
                
                if ( ~isfield(cinfo, 'PixelSpacing') )
                    cinfo.PixelSpacing(1) = 0.779296875;
                    cinfo.PixelSpacing(2) = 0.779296875;
                    %continue;
                end
                
                if ( ~isfield(cinfo, 'SliceThickness') )
                    cinfo.SliceThickness = 2;
                    %continue;
                end
                
                if ( isfield(cinfo, 'InversionTime') )
                    info.InversionTime = cinfo.InversionTime;                
                end
                
                info.ImageSize = [cinfo.Columns cinfo.Rows]; % width, height
                info.PixelSpacing(1) = cinfo.PixelSpacing(1);
                info.PixelSpacing(2) = cinfo.PixelSpacing(2);
                %info.PixelSpacing(3) = cinfo.SpacingBetweenSlices;
                info.PixelSpacing(3) = cinfo.SliceThickness;
                fileInfo = [fileInfo; info];
            end
        end
    end
end

numFile = numel(fileInfo);

disp(['Total ' num2str(numFile) ' is found ...'])

disp('Sorting by patient name ... ')

PatientNames = cell(numFile, 1);
for k=1:numFile
    PatientNames{k} = fileInfo(k).PatientName;
end

[sPatientNames, index] = sort(PatientNames);
fileInfo = adjustOrderInfo(fileInfo, index);

allPatients = cell(0);
allPatients = sPatientNames(1);

for k=2:numFile
    index = [];
    numFile2 = numel(allPatients);
    for p=1:numFile2
        if ( strcmp(PatientNames{k}, allPatients{p}) )
            index = [index p];
        end
    end
    
    if ( isempty(index) )
       allPatients = [allPatients; PatientNames{k} ]; 
    end
end

allPatients{:}

numPatients = numel(allPatients);

headers = [];

for pp = 1:numPatients
    
    patientName = allPatients{pp} % patientName
    index = [];
    numFile = numel(sPatientNames);
    for k=1:numFile
        if ( strcmp(patientName, sPatientNames{k}) )
            index = [index k];
        end
    end
    
    disp('Sorting by serial name ... ')
    cFileInfo = fileInfo(index(:));
    numFile = numel(cFileInfo);
    serialNames = zeros(numFile, 1);
    for k=1:numFile
        if ( ~isempty(cFileInfo(k).SeriesNumber) )
            serialNames(k) = cFileInfo(k).SeriesNumber;
        else
            serialNames(k) = floor((k-1)/10);
        end
    end
    
    [cSerialNames, indexSerial] = sort(serialNames);
    cFileInfo = adjustOrderInfo(cFileInfo, indexSerial);

    allDiffSerials = findDiffItem(cSerialNames);
    
    disp('imaging serials are ')
    allDiffSerials
    
    numOfSerials = numel(allDiffSerials);
    
    for serials = 1:numOfSerials
        
        cSerials = allDiffSerials(serials); % current serials
        
        index = [];
        numFile = numel(cFileInfo);
        for k=1:numFile
            if ( ~isempty(cFileInfo(k).SeriesNumber) )
                seriesNumber = cFileInfo(k).SeriesNumber;
            else
                seriesNumber = floor((k-1)/10);
            end
            
            if ( cSerials == seriesNumber )
                index = [index k];
            end
        end
        cFileInfo_Serials = cFileInfo(index(:));
        
        disp('Sorting by slice position name ... ')
        numFile = numel(cFileInfo_Serials);
        SliceLocation = zeros(numFile, 1);
        for k=1:numFile
            SliceLocation(k) = cFileInfo_Serials(k).SliceLocation;
        end

        [cSliceLocation, indexSlicePosition] = sort(SliceLocation);
        cFileInfo_Serials = adjustOrderInfo(cFileInfo_Serials, indexSlicePosition);

        allDiffSlicePositions = findDiffItem(cSliceLocation);
        allDiffSlicePositions
        
        numOfSlicePositions = numel(allDiffSlicePositions);
        
        for sliceP = 1:numOfSlicePositions
            
            cSlicePosition = allDiffSlicePositions(sliceP);
            
            index = [];
            numFile = numel(cFileInfo_Serials);
            for k=1:numFile
                if ( cSlicePosition == cFileInfo_Serials(k).SliceLocation )
                    index = [index k];
                end
            end
            cFileInfo_SlicePostions = cFileInfo_Serials(index(:));
        
            disp('Sorting by image num name ... ')
            numFile = numel(cFileInfo_SlicePostions);
            if ( numFile < 1 )
                continue;
            end
            
            AcquisitionNumber = zeros(numFile, 1);
            for k=1:numFile
                AcquisitionNumber(k) = cFileInfo_SlicePostions(k).AcquisitionNumber;
            end

            [cAcquisitionNumber, indexAcquisitionNumber] = sort(AcquisitionNumber);
            cFileInfo_SlicePostions = adjustOrderInfo(cFileInfo_SlicePostions, indexAcquisitionNumber);
            
            % save as analyze
            patientName(find(patientName==' ')) = '_';
            patientName(find(patientName=='*')) = 'p';
            Prefix = [patientName '-' num2str(cSerials) '-' num2str(cSlicePosition) '-' num2str(numFile)];
            filename = fullfile(dstDir, [Prefix '.hdr']);
            
            width = cFileInfo_SlicePostions(1).ImageSize(1);
            height = cFileInfo_SlicePostions(1).ImageSize(2);
            depth = numFile;
            
            xvoxelsize = cFileInfo_SlicePostions(1).PixelSpacing(1);
            yvoxelsize = cFileInfo_SlicePostions(1).PixelSpacing(2);
            zvoxelsize = cFileInfo_SlicePostions(1).PixelSpacing(3);
            
%             header = struct('xsize', width, 'ysize', height, 'zsize', depth, ...
%                 'spacingX', xvoxelsize, 'spacingY', yvoxelsize, 'spacingZ', zvoxelsize, ...
%                 'bytes', 2, 'TI', [], 'fileName', [], 'SliceLocation', []);
            
            header = struct('sizeX', width, 'sizeY', height, 'sizeZ', depth, 'sizeT', 1, 'sizeN', 1, 'sizeM', 1, ... 
                            'spacingX', xvoxelsize, 'spacingY', yvoxelsize, 'spacingZ', zvoxelsize, ... 
                            'spacingT', 1.0, 'spacingN', 1.0, 'spacingM', 1.0, 'positionPatient', [0 0 0], 'orientationPatient', eye(3,3));
    
            header.TI = zeros(numFile, 1);
            header.SliceLocation = zeros(numFile, 1);
            header.fileName = cell(numFile, 1);
            
            data = zeros([height width depth], 'uint32');
            
            for k=1:numFile
                disp(['reading image ' cFileInfo_SlicePostions(k).filename]);
                sliceData = dicomread(cFileInfo_SlicePostions(k).filename);
                if ( size(sliceData, 3) > 0 )
                    data(:, :, k) = uint32(mean(double(sliceData), 3));
                else
                    data(:, :, k) = sliceData;
                end
                
                header.TI(k) = cFileInfo_SlicePostions(k).InversionTime;
                header.fileName{k} = cFileInfo_SlicePostions(k).filename;
                header.SliceLocation(k) = cFileInfo_SlicePostions(k).SliceLocation;
            end
            header2 = CreateFtkHeaderInfoFrom(data, header);
            Matlab_SaveAnalyze(single(data), header2, filename);
            
            headers = [headers header];
        end
    end
end


