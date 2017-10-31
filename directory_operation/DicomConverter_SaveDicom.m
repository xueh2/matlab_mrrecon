% function DicomConverter_SaveDicom(home, dstDir)

% home = 'D:\data\4DRecon_20111206\newVolunteers\Dicom_v4\Dicom';
% dstDir = 'D:\data\4DRecon_20111206\newVolunteers\Volunteer4';

% home = 'D:\data\AnwarWB_follow_up\composed\CL\dicom';
% dstDir = 'D:\data\AnwarWB_follow_up\composed\CL\SeperateSeries';
 
home = 'D:\data\AnwarWB_follow_up\composed\SA\dicom';
dstDir = 'D:\data\AnwarWB_follow_up\composed\SA\SeperateSeries';

disp('parsing the home directory ...')

subdirs = [];
subdirs = FindDeepSubDirs(home, subdirs, [])

% subdirs = [subdirs {home}];

info = struct('filename', cell(0), 'PatientName', '', 'SeriesNumber', -1, 'SliceLocation', -1,...
    'AcquisitionNumber', -1, 'ImageSize', [-1 -1], 'PixelSpacing', [0 0 0]);
fileInfo = [];

SeriesNumbers = [];

num = numel(subdirs);
for k=1:num
    
    D = dir( fullfile( subdirs{k} ,'*.*') );
    fileNum = numel(D);
    
    disp(['Read images at ' subdirs{k}])
    
    for s = 1:fileNum
        
        if ( mod(s, 1000) == 0 )
            disp([num2str(s) ' ... ']);
        end
        
        if ( D(s).isdir == 0 )
            
            currentFileName = fullfile(subdirs{k}, D(s).name);
            if ( isdicom(currentFileName) )
%                 info = struct('filename', '', 'PatientName', '', 'SeriesNumber', -1, 'SliceLocation', -1,...
%                     'AcquisitionNumber', -1, 'ImageSize', [-1 -1], 'PixelSpacing', [0 0 0]);

                cinfo = dicominfo(currentFileName);
%                 info.filename = currentFileName;
                
%                 if ( isfield(cinfo, 'PatientName') )
%                     info.PatientName = cinfo.PatientName.FamilyName;
%                 else
%                     info.PatientName = 'Unknowun';
%                 end
                
                if ( isfield(cinfo, 'SeriesNumber') )
                    % info.SeriesNumber = cinfo.SeriesNumber;
                    ind = find(SeriesNumbers == cinfo.SeriesNumber);
                    if ( sum(ind) == 0 )
                        SeriesNumbers = [SeriesNumbers cinfo.SeriesNumber];
                    end
                end
                                
                
%                 if ( isfield(cinfo, 'SliceLocation') )
%                     info.SliceLocation = cinfo.SliceLocation;
%                 end
%                 
%                 if ( isfield(cinfo, 'AcquisitionNumber') )
%                     info.AcquisitionNumber = cinfo.AcquisitionNumber;                
%                 end
% 
%                 if ( isfield(cinfo, 'ImageNumber') )
%                     info.AcquisitionNumber = cinfo.ImageNumber;                
%                 end
%                 
%                 if ( ~isfield(cinfo, 'PixelSpacing') )
%                     cinfo.PixelSpacing(1) = 0.779296875;
%                     cinfo.PixelSpacing(2) = 0.779296875;
%                     %continue;
%                 end
%                 
%                 if ( ~isfield(cinfo, 'SliceThickness') )
%                     cinfo.SliceThickness = 2;
%                     %continue;
%                 end
%                 
%                 info.ImageSize = [cinfo.Columns cinfo.Rows]; % width, height
%                 info.PixelSpacing(1) = cinfo.PixelSpacing(1);
%                 info.PixelSpacing(2) = cinfo.PixelSpacing(2);
%                 %info.PixelSpacing(3) = cinfo.SpacingBetweenSlices;
%                 info.PixelSpacing(3) = cinfo.SliceThickness;
%                 fileInfo = [fileInfo; info];
            end
        end
    end
end

numSeries = numel(SeriesNumbers);
folders = cell(numSeries, 1);
for i=1:numSeries
    folders{i} = fullfile(dstDir, ['Series' num2str(SeriesNumbers(i))]);
    mkdir(folders{i});
end

num = numel(subdirs);
for k=1:num
    
    D = dir( fullfile( subdirs{k} ,'*.*') );
    fileNum = numel(D);
    
    disp(['Read images at ' subdirs{k}])
    
    for s = 1:fileNum
        
        if ( mod(s, 1000) == 0 )
            disp([num2str(s) ' ... ']);
        end
        
        if ( D(s).isdir == 0 )
            
            currentFileName = fullfile(subdirs{k}, D(s).name);
            if ( isdicom(currentFileName) )
                cinfo = dicominfo(currentFileName);
                
                if ( isfield(cinfo, 'SeriesNumber') )
                    ind = find(SeriesNumbers==cinfo.SeriesNumber);
                    copyfile(currentFileName, fullfile(folders{ind}, D(s).name));
                end
                
            end
        end
    end
end

% disp('Sorting by patient name ... ')
% 
% PatientNames = cell(numFile, 1);
% for k=1:numFile
%     PatientNames{k} = fileInfo(k).PatientName;
% end
% 
% [sPatientNames, index] = sort(PatientNames);
% fileInfo = adjustOrderInfo(fileInfo, index);
% 
% allPatients = cell(0);
% allPatients = sPatientNames(1);
% 
% for k=2:numFile
%     index = [];
%     numFile2 = numel(allPatients);
%     for p=1:numFile2
%         if ( strcmp(PatientNames{k}, allPatients{p}) )
%             index = [index p];
%         end
%     end
%     
%     if ( isempty(index) )
%        allPatients = [allPatients; PatientNames{k} ]; 
%     end
% end
% 
% allPatients{:}
% 
% numPatients = numel(allPatients);
% 
% for pp = 1:numPatients
%     
%     patientName = allPatients{pp} % patientName
%     index = [];
%     numFile = numel(sPatientNames);
%     for k=1:numFile
%         if ( strcmp(patientName, sPatientNames{k}) )
%             index = [index k];
%         end
%     end
%     
%     disp('Sorting by serial name ... ')
%     cFileInfo = fileInfo(index(:));
%     numFile = numel(cFileInfo);
%     serialNames = zeros(numFile, 1);
%     for k=1:numFile
%         if ( ~isempty(cFileInfo(k).SeriesNumber) )
%             serialNames(k) = cFileInfo(k).SeriesNumber;
%         else
%             serialNames(k) = floor((k-1)/10);
%         end
%     end
%     
%     [cSerialNames, indexSerial] = sort(serialNames);
%     cFileInfo = adjustOrderInfo(cFileInfo, indexSerial);
% 
%     allDiffSerials = findDiffItem(cSerialNames);
%     
%     disp('imaging serials are ')
%     allDiffSerials
%     
%     numOfSerials = numel(allDiffSerials);
%     
%     for serials = 1:numOfSerials
%         
%         cSerials = allDiffSerials(serials); % current serials
%         
%         index = [];
%         numFile = numel(cFileInfo);
%         for k=1:numFile
%             if ( ~isempty(cFileInfo(k).SeriesNumber) )
%                 seriesNumber = cFileInfo(k).SeriesNumber;
%             else
%                 seriesNumber = floor((k-1)/10);
%             end
%             
%             if ( cSerials == seriesNumber )
%                 index = [index k];
%             end
%         end
%         cFileInfo_Serials = cFileInfo(index(:));
%         
%         disp('Sorting by slice position name ... ')
%         numFile = numel(cFileInfo_Serials);
%         SliceLocation = zeros(numFile, 1);
%         for k=1:numFile
%             SliceLocation(k) = cFileInfo_Serials(k).SliceLocation;
%         end
% 
%         [cSliceLocation, indexSlicePosition] = sort(SliceLocation);
%         cFileInfo_Serials = adjustOrderInfo(cFileInfo_Serials, indexSlicePosition);
% 
%         allDiffSlicePositions = findDiffItem(cSliceLocation);
%         allDiffSlicePositions
%         
%         numOfSlicePositions = numel(allDiffSlicePositions);
%         
%         for sliceP = 1:numOfSlicePositions
%             
%             cSlicePosition = allDiffSlicePositions(sliceP);
%             
%             index = [];
%             numFile = numel(cFileInfo_Serials);
%             for k=1:numFile
%                 if ( cSlicePosition == cFileInfo_Serials(k).SliceLocation )
%                     index = [index k];
%                 end
%             end
%             cFileInfo_SlicePostions = cFileInfo_Serials(index(:));
%         
%             disp('Sorting by image num name ... ')
%             numFile = numel(cFileInfo_SlicePostions);
%             if ( numFile < 1 )
%                 continue;
%             end
%             
%             AcquisitionNumber = zeros(numFile, 1);
%             for k=1:numFile
%                 AcquisitionNumber(k) = cFileInfo_SlicePostions(k).AcquisitionNumber;
%             end
% 
%             [cAcquisitionNumber, indexAcquisitionNumber] = sort(AcquisitionNumber);
%             cFileInfo_SlicePostions = adjustOrderInfo(cFileInfo_SlicePostions, indexAcquisitionNumber);
%             
%             % save as analyze
%             patientName(find(patientName==' ')) = '_';
%             patientName(find(patientName=='*')) = 'p';
%             Prefix = [patientName '-' num2str(cSerials) '-' num2str(cSlicePosition) '-' num2str(numFile)];
%             filename = fullfile(dstDir, [Prefix '.hdr']);
%             
%             width = cFileInfo_SlicePostions(1).ImageSize(1);
%             height = cFileInfo_SlicePostions(1).ImageSize(2);
%             depth = numFile;
%             
%             xvoxelsize = cFileInfo_SlicePostions(1).PixelSpacing(1);
%             yvoxelsize = cFileInfo_SlicePostions(1).PixelSpacing(2);
%             zvoxelsize = cFileInfo_SlicePostions(1).PixelSpacing(3);
%             
%             header = struct('xsize', width, 'ysize', height, 'zsize', depth, ...
%                 'xvoxelsize', xvoxelsize, 'yvoxelsize', yvoxelsize, 'zvoxelsize', zvoxelsize, ...
%                 'bytes', 2);
%             data = zeros([height width depth], 'uint32');
%             
%             for k=1:numFile
%                 disp(['reading image ' cFileInfo_SlicePostions(k).filename]);
%                 sliceData = dicomread(cFileInfo_SlicePostions(k).filename);
%                 if ( size(sliceData, 3) > 0 )
%                     data(:, :, k) = uint32(mean(double(sliceData), 3));
%                 else
%                     data(:, :, k) = sliceData;
%                 end
%             end
%             header = CreateFtkHeaderInfoFrom(data, header);
%             Matlab_SaveAnalyze(data, header, filename);
%         end
%     end
% end


