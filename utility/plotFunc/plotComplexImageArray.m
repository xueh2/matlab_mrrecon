
function [mag, header] = plotComplexImageArray(im, voxelsize, newNpe, centre, width, delayTime, endo_pt, epi_pt, plotFlag)

if ( nargin < 2 | isempty(voxelsize) )
    voxelsize = [1 1 1];
end

if ( nargin < 3 | isempty(newNpe) )
    newNpe = size(im,2);
end

if ( nargin < 5 )
    centre = 1024;
    width = 1024;
end

if ( nargin < 6 )
   delayTime = 1/30;
end

if( nargin < 7 )
    endo_pt = [];
end

if( nargin < 8 )
    epi_pt = [];
end

if ( nargin < 9 )
   plotFlag = 1;
end

s = size(im);

smallVoxelSize = min(voxelsize(1), voxelsize(2));

if ( length(s) == 4 )
    Nfe = size(im, 1);
    Npe = size(im, 2);
    numOfCoils = size(im, 3);
    numOfFrames = size(im, 4);

    Img = SoS_Image_TemporalArray(im);

    if ( newNpe > Npe )
        Img = Zero_Padding_Resize_NoFiltering(Img, Nfe, newNpe);
    end
    
    %Img = flipdim(Img, 2);    
    header = CreateFtkHeaderInfo(Img, voxelsize);

    mag = abs(Img); 
    lowR = 0;
    highR = 4096;
    mag = normalizeImage2Range(mag, lowR, highR);
        
%     % shift by 1 pixel to match the inline icepat images
%     A = [1 0 1; 0 1 0; 0 0 1];
%     tform = maketform('affine',A');
%     for f=1:numOfFrames
%         s = imtransform(mag(:,:,f),tform, 'nearest', 'XData',[1 size(mag,2)],'YData',[1 size(mag,1)]);
%         mag(:,:,f) = s;
%     end
    
    if ( plotFlag )
        magWindowed = normalizeWindowSetting(mag, centre, width);        
        imArrayMrFtkPlayer(magWindowed, header,delayTime);
    end
end

if ( length(s) == 3 )
    Nfe = size(im, 1);
    Npe = size(im, 2);
    numOfFrames = size(im, 3);

    if ( newNpe > Npe )
        im = Zero_Padding_Resize_NoFiltering(im, Nfe, newNpe);
    end
    
    header = CreateFtkHeaderInfo(im, voxelsize);
    %im = flipdim(im, 2);

    mag = abs(im); 
    lowR = 0;
    highR = 4096;
    mag = normalizeImage2Range(mag, lowR, highR);

%     % shift by 1 pixel to match the inline icepat images
%     A = [1 0 1; 0 1 0; 0 0 1];
%     tform = maketform('affine',A');
%     for f=1:numOfFrames
%         s = imtransform(mag(:,:,f),tform, 'nearest', 'XData',[1 size(mag,2)],'YData',[1 size(mag,1)]);
%         mag(:,:,f) = s;
%     end

    if ( plotFlag )
        magWindowed = normalizeWindowSetting(mag, centre, width);        
        imArrayMrFtkPlayer(magWindowed, header,delayTime, endo_pt, epi_pt);
    end
end

if ( length(s) == 2 )
    Nfe = size(im, 1);
    Npe = size(im, 2);

    if ( newNpe > Npe )
        im = Zero_Padding_Resize_NoFiltering(im, Nfe, newNpe);
    end
    
    header = CreateFtkHeaderInfo(im, voxelsize);
    %im = flipdim(im, 2);

    mag = abs(im); 
    lowR = 0;
    highR = 4096;
    mag = normalizeImage2Range(mag, lowR, highR);

%     % shift by 1 pixel to match the inline icepat images
%     A = [1 0 1; 0 1 0; 0 0 1];
%     tform = maketform('affine',A');
%     for f=1:numOfFrames
%         s = imtransform(mag(:,:,f),tform, 'nearest', 'XData',[1 size(mag,2)],'YData',[1 size(mag,1)]);
%         mag(:,:,f) = s;
%     end

    if ( plotFlag )
        magWindowed = normalizeWindowSetting(mag, centre, width);        
        imArrayMrFtkPlayer(magWindowed, header,delayTime, endo_pt, epi_pt);
    end
end

header.spacingX = smallVoxelSize;
header.spacingY = smallVoxelSize;
