
function [mag, magResized, header] = plotComplexImageArrayMontage(im, voxelsize, centre, width, sizeRatio, nRows)


if ( nargin < 2 )
    voxelsize = [1 1 1];
end

if ( nargin < 4 )
    centre = 1024;
    width = 1024;
end

if ( nargin < 5 )
    sizeRatio = 2;
end

if ( nargin < 6 )
    nRows = 1;
end

s = size(im);

smallVoxelSize = min(voxelsize(1), voxelsize(2));

if ( length(s) == 4 )
    Nfe = size(im, 1);
    Npe = size(im, 2);
    numOfCoils = size(im, 3);
    numOfFrames = size(im, 4);

    Img = SoS_Image_TemporalArray(im);

    header = CreateFtkHeaderInfo(Img, voxelsize);

    lowR = 0;
    highR = 4096;
    mag = abs(Img); 
    mag = normalizeImage2Range(mag, lowR, highR);
    mag = normalizeWindowSetting(mag, centre, width);
    magResized = resizeImageVolume(mag, header, sizeRatio);
    % imArrayplayer(magResized,delayTime);
    s = size(magResized);
    D = reshape(magResized, [s(1) s(2) 1 s(3)]);
    figure;montage(uint8(D), [nRows, NaN]); 
end

if ( length(s) == 3 )
    Nfe = size(im, 1);
    Npe = size(im, 2);
    numOfFrames = size(im, 3);

    header = CreateFtkHeaderInfo(im, voxelsize);

    lowR = 0;
    highR = 4096;
    mag = abs(im); 
    mag = normalizeImage2Range(mag, lowR, highR);
    mag = normalizeWindowSetting(mag, centre, width);
    magResized = resizeImageVolume(mag, header, sizeRatio);
    %imArrayplayer(magResized,delayTime);
    s = size(magResized);
    D = reshape(magResized, [s(1) s(2) 1 s(3)]);
    figure;montage(uint8(D), 'Size', [nRows, NaN]);     
end

header.spacingX = smallVoxelSize;
header.spacingY = smallVoxelSize;
