function DicomFixTriggerTime(home, dstDir, TriggerTime)

mkdir(dstDir);

D = dir( fullfile( home,'*.*') );
fileNum = numel(D);

disp(['Read images at ' home])
    
for s = 1:fileNum
    if ( D(s).isdir == 0 )

        currentFileName = fullfile(home, D(s).name)
        if ( isdicom(currentFileName) )
            info = struct('filename', '', 'PatientName', '', 'SeriesNumber', -1, 'SliceLocation', -1,...
                'AcquisitionNumber', -1, 'AcquisitionTime', -1, 'ImageSize', [-1 -1], 'PixelSpacing', [0 0 0]);

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

            if ( isfield(cinfo, 'AcquisitionTime') )
                info.AcquisitionTime = cinfo.AcquisitionTime;                
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

            info.ImageSize = [cinfo.Columns cinfo.Rows]; % width, height
            info.PixelSpacing(1) = cinfo.PixelSpacing(1);
            info.PixelSpacing(2) = cinfo.PixelSpacing(2);
            %info.PixelSpacing(3) = cinfo.SpacingBetweenSlices;
            info.PixelSpacing(3) = cinfo.SliceThickness;
            
            tt = TriggerTime(cinfo.InstanceNumber);            
            cinfo.TriggerTime = tt;            
            cinfo.SOPInstanceUID = sprintf('%s%s',cinfo.SOPInstanceUID ,'1');
            cinfo.SeriesInstanceUID = sprintf('%s%s',cinfo.SeriesInstanceUID,'1');
            cinfo.SeriesNumber = cinfo.SeriesNumber * 10;
            
            filename = fullfile(dstDir, [num2str(cinfo.InstanceNumber) '.IMA']);
            
            d = dicomread(currentFileName);            
            dicomwrite(int16(d), filename, cinfo);            
        end
    end
end
