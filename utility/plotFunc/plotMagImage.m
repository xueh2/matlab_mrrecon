
function plotMagImage(im, voxelsize, centre, width, sizeRatio, delayTime)

Nfe = size(im, 1);
Npe = size(im, 2);
numOfCoils = size(im, 3);

header = CreateFtkHeaderInfo(im, voxelsize);

lowR = 0;
highR = 4096;
mag = im; 
mag = normalizeImage2Range(mag, lowR, highR);
mag = normalizeWindowSetting(mag, centre, width);
% magResized = resizeImageVolume(mag, header, sizeRatio);
plotMrFtkImage(mag, header, -1, 1);