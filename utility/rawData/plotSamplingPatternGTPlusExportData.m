
function plotSamplingPatternGTPlusExportData(data)

RO = size(data, 1);
E1 = size(data, 2);
CHA = size(data, 3);
E2 = size(data, 5);

NDim = numel(size(data));

mag = abs(data);
for ii=NDim:-1:6
    mag = sum(mag, ii);
end

mag2 = sum(mag, 5);
plotKSpaceSamplingPattern(squeeze(mag2(:,:,1, 1)));
title('RO-E1 plane');

mag2 = sum(mag, 1);
if ( E2 > 1)
    plotKSpaceSamplingPattern(squeeze(mag2(1,:,1,1,:)));
    title('E1-E2 plane');
end
