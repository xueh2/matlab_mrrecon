
function plotComplexImage(im, voxelsize, centre, width)

if ( nargin < 2 )
    voxelsize = [1 1 1];
end

if ( nargin < 4 )
    centre = 250;
    width = 500;
end

Nfe = size(im, 1);
Npe = size(im, 2);
numOfCoils = size(im, 3);

Img = SoS_Image(im);

header = CreateFtkHeaderInfo(Img, voxelsize);

lowR = 0;
highR = 4096;
mag = abs(Img); 
mag = normalizeImage2Range(mag, lowR, highR);
mag = normalizeWindowSetting(mag, centre, width);
% magResized = resizeImageVolume(mag, header, sizeRatio);
% figure;imshow(abs(magResized), [], 'Border', 'tight');
% imtool(abs(magResized), []);
plotMrFtkImage(mag, header, -1, 1);