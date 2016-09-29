function [keyFrame, mssd] = PerformKeyFrameSelectionSSD(data)
% perform motion correction using SSD to pick a key frame

numOfIm = size(data, 3);
ssd = zeros(numOfIm, numOfIm);
for ff=1:numOfIm
    for ff2=ff:numOfIm
        diffIm = data(:,:,ff)-data(:,:,ff2);
        % ssd(ff, ff2) = norm(diffIm(:));
        diffIm = diffIm .* diffIm;
        ssd(ff, ff2) = sqrt(sum(diffIm(:)));
        if ( ff ~= ff2 ) 
            ssd(ff2, ff) = ssd(ff, ff2);
        end
    end
end

ssd = sort(ssd, 1);
mssd = median(ssd);
[minssd, keyFrame] = min(mssd);
keyFrame = keyFrame - 1;
