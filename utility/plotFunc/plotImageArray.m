
function [mag, header] = plotImageArray(im, voxelsize, centre, width, delayTime, plotFlag)

if ( nargin < 2 )
    voxelsize = [1 1 1];
end

if ( nargin < 4 )
    centre = 1024;
    width = 1024;
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
    
    header = CreateFtkHeaderInfo(im, voxelsize);

    lowR = 0;
    highR = 4096;
    mag = abs(im); 
    mag = normalizeImage2Range(mag, lowR, highR);

    if ( plotFlag )
        magWindowed = normalizeWindowSetting(mag, centre, width);        
        imArrayMrFtkPlayer(magWindowed, header,delayTime);
    end
end

header.spacingX = smallVoxelSize;
header.spacingY = smallVoxelSize;
