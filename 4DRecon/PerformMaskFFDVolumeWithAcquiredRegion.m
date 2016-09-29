function [volumeMasked, headerMasked, mask] = PerformMaskFFDVolumeWithAcquiredRegion(volume, header, imageSize, imagePositionPatient, imageOrientationPatient, linearOrRotate)
% mask the FFD interpolated volume with acquired region; if not in the acquired region, the pixel is set to be zero
% volume : [COL LIN SLC PHS], reconed 4D volume
% imageSize : image size for original 2D frame
% image position and orientation for every 2D frame
% the order the 2D frames are sorted
% linearOrRotate : 'linear' for linear sweeping; 'rotate' for rotating acquisition

volumeMasked = volume;
headerMasked = header;
mask = zeros(size(volume));

if ( strcmp(linearOrRotate, 'linear') == 1 )
    
    % composite a header for the acquired region
    headerAcquired = header;
    headerAcquired.sizeX = imageSize(1);
    headerAcquired.sizeY = imageSize(2);
    headerAcquired.sizeZ = size(imagePositionPatient, 1);
    
    distZ = norm(imagePositionPatient(end, :)-imagePositionPatient(1, :));
    headerAcquired.spacingZ = distZ/headerAcquired.sizeZ;
    
    headerAcquired.positionPatient = imagePositionPatient(1, :);
    headerAcquired.orientationPatient(1, :) = imageOrientationPatient(1, 1:3);
    headerAcquired.orientationPatient(2, :) = imageOrientationPatient(1, 4:6);
    headerAcquired.orientationPatient(3, :) = cross(headerAcquired.orientationPatient(1,:), headerAcquired.orientationPatient(2,:));
    
    for z=1:header.sizeZ
        z
        for y=1:header.sizeY
            for x=1:header.sizeX                
                [wx, wy, wz ] = Image2WorldMrFtk(header, x-1, y-1, z-1);
                inRange = IsInRangeMrFtk(headerAcquired, wx, wy, wz);
                if ( ~inRange )
                    volumeMasked(y, x, z) = 0;
                else
                    mask(y, x, z) = 1;
                end
            end
        end
    end    
end

if ( strcmp(linearOrRotate, 'rotate') == 1 )
    
    % compute the rotating axis   
    N = size(imagePositionPatient, 1);
    axisStart = zeros(N, 3);
    axisEnd = zeros(N, 3);
    radius = zeros(N, 1);
    
    for n=1:N        
        header2D = Dicom2HeaderMrFtk(zeros(imageSize(1), imageSize(2)), [header.spacingX header.spacingY header.spacingZ], imagePositionPatient(n, :), imageOrientationPatient(n, :));
        [wc1, wc2, wc3, wc4] = computeFourCorner(zeros(imageSize(1), imageSize(2)), header2D);        
        axisStart(n, :) = (wc1+wc4)/2;
        axisEnd(n, :) = (wc2+wc3)/2;
        radius(n) = norm(wc1-wc4)/2;
    end
    
    figure;
    hold on
    plot3(imagePositionPatient(:,1), imagePositionPatient(:,2), imagePositionPatient(:,3), '+')
    plot3(axisStart(:,1), axisStart(:,2), axisStart(:,3), 'r+');
    plot3(axisEnd(:,1), axisEnd(:,2), axisEnd(:,3), 'k+');
    hold off
    
    sPt = mean(axisStart);
    ePt = mean(axisEnd);    
    r = mean(radius);
    
    A = ePt(1) - sPt(1);
    B = ePt(2) - sPt(2);
    C = ePt(3) - sPt(3);
    normL = (A*A+B*B+C*C);
    
    for z=1:header.sizeZ
        z
        for y=1:header.sizeY
            for x=1:header.sizeX                
    
                [wx, wy, wz ] = Image2WorldMrFtk(header, x-1, y-1, z-1);
                
                pt = [wx wy wz];
                t = A*(sPt(1)-wx)+B*(sPt(2)-wy)+C*(sPt(3)-wz);
                t = -t/normL;
                
                crossPt = [sPt(1)+A*t sPt(2)+B*t sPt(3)+C*t];
                
                if ( t<0 | t>1 | norm(crossPt-pt)>r )
                    volumeMasked(y, x, z) = 0;
                else
                    mask(y, x, z) = 1;
                end
            end
        end
    end
end
