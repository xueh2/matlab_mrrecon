function Render4DVolume(binData, pixelSpacing, ImagePositionPatient, ImageOrientationPatient, slc, phs, plotStep, centre, width, delayTime, saveAvi, aviName)
% Render the 4D volume for selected slice or phase
% if slc > 0, render slice slc; if phs > 0, render phase phs

if ( ~exist('centre') )
    centre = 1024;
end

if ( ~exist('width') )
    width = 1024;
end

if ( ~exist('delayTime') )
    delayTime = 1/20;
end

if ( ~exist('saveAvi') )
    saveAvi = 0;
end

if ( ~exist('aviName') )
    aviName = 'd:/temp/3D.avi';
end

if ( ~exist('plotStep') )
    plotStep = 1;
end

if ( slc > 0 )
    
    data = squeeze(binData(:, :, slc, :));
    
    lowR = 0;
    highR = 4096;
    mag = normalizeImage2Range(data, lowR, highR);
    magWindowed = normalizeWindowSetting(mag, centre, width);        
        
    header = Dicom2HeaderMrFtk(magWindowed, pixelSpacing, ImagePositionPatient(1,:), ImageOrientationPatient(1,:));
    
    imArrayMrFtkPlayer(magWindowed, header, delayTime);    
end

if ( phs > 0 )
   
    data = squeeze(binData(:, :, :, phs));   
    num = size(data, 3);
    
    if ( saveAvi )
        vidObj = VideoWriter(aviName, 'Motion JPEG AVI');
        open(vidObj);
        ha = -1;
        for k=1:plotStep:num
            d = data(:,:,k);
            h = Dicom2HeaderMrFtk(d, pixelSpacing, ImagePositionPatient(k,:), ImageOrientationPatient(k,:));
            ha = plotMrFtkImageIn3D(d, h, ha, centre, width);
            frame = getframe(gcf);
            writeVideo(vidObj,frame);
            % aviobj = addframe(aviobj,frame);
        end
        close(vidObj);
        %aviobj = close(aviobj);
    else        
        ha = -1;
        for k=1:plotStep:num
            d = data(:,:,k);
            h = Dicom2HeaderMrFtk(d, pixelSpacing, ImagePositionPatient(k,:), ImageOrientationPatient(k,:));
            ha = plotMrFtkImageIn3D(d, h, ha, centre, width);
        end
    end
end