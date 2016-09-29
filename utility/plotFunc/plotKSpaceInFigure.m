
function h = plotKSpaceInFigure(kspace, voxelsize, centre, width, sizeRatio, delayTime, h)

Nfe = size(kspace, 1);
Npe = size(kspace, 2);
numOfCoils = size(kspace, 3);

Img = SoS(kspace);

header = CreateFtkHeaderInfo(Img, voxelsize);

lowR = 0;
highR = 4096;
mag = abs(Img); 
mag = normalizeImage2Range(mag, lowR, highR);
mag = normalizeWindowSetting(mag, centre, width);
% magResized = resizeImageVolume(mag, header, sizeRatio);
% figure;imshow(abs(magResized), []);
% imtool(abs(magResized), []);
[h, imH] = plotMrFtkImage(mag, header, h, 1);