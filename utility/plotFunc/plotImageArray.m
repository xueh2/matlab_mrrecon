
function [mag, header, xyAxes] = plotImageArray(im, voxelsize, landmarks, centre, width, delayTime, plotFlag)
% [mag, header, xyAxes] = plotImageArray(im, voxelsize, landmarks, centre, width, delayTime, plotFlag)

if ( nargin < 2 )
    voxelsize = [1 1 1];
end

if ( nargin < 3 )
    landmarks = [];
end

if ( nargin < 5 )
    centre = 1024;
    width = 1024;
end

if ( nargin < 6 )
   delayTime = 0.1;
end

if ( nargin < 7 )
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
        xyAxes = imArrayMrFtkPlayer(magWindowed, header,delayTime, [], [], landmarks);
    end
end

header.spacingX = smallVoxelSize;
header.spacingY = smallVoxelSize;
