function DicomConverter(home)

disp('parsing the home directory ...')

subdirs = [];
subdirs = FindDeepSubDirs(home, subdirs, [])

% subdirs = [subdirs {home}];

info = struct('filename', cell(0), 'PatientName', '', 'SeriesNumber', -1, 'SliceLocation', -1,...
    'AcquisitionNumber', -1, 'ImageSize', [-1 -1 -1], 'PixelSpacing', [0 0 0]);
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
                    'AcquisitionNumber', -1, 'ImageSize', [-1 -1], 'PixelSpacing', [0 0 0]);

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
                
                info.ImageSize = [cinfo.Columns cinfo.Rows]; % width, height
                info.PixelSpacing = cinfo.PixelSpacing;
                
                fileInfo = [fileInfo; info];
            end
        end
    end
end

numFile = numel(fileInfo);

disp(['Total ' num2str(numFile) ' is found ...'])

disp('Sorting by patient name ... ')

disp('Sorting by study name ... ')

disp('Sorting by slice position name ... ')

disp('Sorting by image num name ... ')
