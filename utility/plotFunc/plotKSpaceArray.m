
function [mag, header] = plotKSpaceArray(kspace, voxelsize, centre, width, delayTime);

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

    Nfe = size(kspace, 1);
    Npe = size(kspace, 2);
    numOfCoils = size(kspace, 3);
    numOfFrames = size(kspace, 4);

    Img = SoS_TemporalArray(kspace);

    header = CreateFtkHeaderInfo(Img, voxelsize);

    lowR = 0;
    highR = 4096;
    mag = abs(Img); 
    mag = normalizeImage2Range(mag, lowR, highR);
    mag = normalizeWindowSetting(mag, centre, width);
    imArrayMrFtkPlayer(mag, header,delayTime);

end

if ( length(s) == 3 )
    Nfe = size(kspace, 1);
    Npe = size(kspace, 2);
    numOfFrames = size(kspace, 3);

    Img = ifft2DForVolume(kspace);
    header = CreateFtkHeaderInfo(Img, voxelsize);

    lowR = 0;
    highR = 4096;
    mag = abs(Img); 
    mag = normalizeImage2Range(mag, lowR, highR);
    mag = normalizeWindowSetting(mag, centre, width);
    imArrayMrFtkPlayer(mag, header,delayTime);
end

if ( length(s) == 5 )

    Nfe = size(kspace, 1);
    Npe = size(kspace, 2);
    numOfCoils = size(kspace, 3);
    numOfFrames = size(kspace, 4);
    numOfSets = size(kspace, 5);
    
    for set=1:s(5)
        Img = SoS_TemporalArray(kspace(:,:,:,:,set));

        header = CreateFtkHeaderInfo(Img, voxelsize);

        lowR = 0;
        highR = 4096;
        mag = abs(Img); 
        mag = normalizeImage2Range(mag, lowR, highR);
        mag = normalizeWindowSetting(mag, centre, width);
        imArrayMrFtkPlayer(mag, header,delayTime);
    end
end
