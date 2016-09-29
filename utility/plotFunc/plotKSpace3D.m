
function [mag, header] = plotKSpace3D(kspace, voxelsize, centre, width, delayTime);
% kspace : COL LIN PAR CHA

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

s = size(kspace);

if ( length(s) == 4 )

    Im = ifft3c(kspace);
    Img = SoS_Image(Im);
    % imagescn(real(Img), [], [], [], 3);

    header = CreateFtkHeaderInfo(Img, voxelsize);

    lowR = 0;
    highR = 4096;
    mag = abs(Img); 
    mag = normalizeImage2Range(mag, lowR, highR);
    mag = normalizeWindowSetting(mag, centre, width);
    imArrayMrFtkPlayer(mag, header,delayTime);

end
