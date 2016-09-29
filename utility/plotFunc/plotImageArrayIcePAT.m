
function [mag, header] = plotImageArrayIcePAT(im, voxelsize, newNpe, centre, width, delayTime, plotFlag)

if ( nargin < 2 )
    voxelsize = [1 1 1];
end

if ( nargin < 3 )
    newNpe = size(im,2);
end

if ( nargin < 5 )
    centre = 200;
    width = 400;
end

if ( nargin < 5 )
   delayTime = 0.1;
end

if ( nargin < 6 )
   plotFlag = 1;
end

s = size(im);

smallVoxelSize = min(voxelsize(1), voxelsize(2));

if ( length(s) == 3 )
    Nfe = size(im, 1);
    Npe = size(im, 2);
    numOfFrames = size(im, 3);

    if ( newNpe > Npe )
        im = Zero_Padding_Resize(im, Nfe, newNpe);
    end
    
    header = CreateFtkHeaderInfo(im, voxelsize);
    im = flipdim(im, 2);

    mag = im;
    
    % shift by 1 pixel to match the inline icepat images
%     A = [1 0 1; 0 1 0; 0 0 1];
%     tform = maketform('affine',A');
%     for f=1:numOfFrames
%         s = imtransform(im(:,:,f),tform, 'nearest', 'XData',[1 size(im,2)],'YData',[1 size(im,1)]);
%         im(:,:,f) = s;
%     end

    if ( plotFlag )
        imWindowed = normalizeWindowSetting(im, centre, width);
        imArrayMrFtkPlayer(imWindowed, header,delayTime);
    end
end

if ( length(s) == 2 )
    Nfe = size(im, 1);
    Npe = size(im, 2);
    numOfFrames = 1;

    if ( newNpe > Npe )
        im = Zero_Padding_Resize(im, Nfe, newNpe);
    end
    
    header = CreateFtkHeaderInfo(im, voxelsize);
    im = flipdim(im, 2);

    mag = im;
    
    % shift by 1 pixel to match the inline icepat images
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

header.spacingX = smallVoxelSize;
header.spacingY = smallVoxelSize;
