
function plotKSpace(kspace, voxelsize, centre, width, sizeRatio)

if ( nargin < 2 )
    voxelsize = [1 1 1];
end

if ( nargin < 4 )
    centre = 1024;
    width = 2048;
end

if ( nargin < 5 )
    sizeRatio = floor(512/size(kspace,1));
    if ( sizeRatio < 1 )
        sizeRatio = 1;
    end
end

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
plotMrFtkImage(mag, header, -1, 1);