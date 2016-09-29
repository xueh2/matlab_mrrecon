
function plotMagImageArray(im, voxelsize, centre, width, sizeRatio, delayTime)

Nfe = size(im, 1);
Npe = size(im, 2);
numOfCoils = size(im, 3);
numOfFrames = size(im, 4);

header = CreateFtkHeaderInfo(im, voxelsize);

lowR = 0;
highR = 4096;
mag = normalizeImage2Range(im, lowR, highR);
mag = normalizeWindowSetting(mag, centre, width);
magResized = resizeImageVolume(mag, header, sizeRatio);
imArrayplayer(magResized,delayTime);