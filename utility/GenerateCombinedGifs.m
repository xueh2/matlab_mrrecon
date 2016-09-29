
function [combinedData, headerCombined] = GenerateCombinedGifs(icepat, vicpaas, header, giffile, sizeRatio, centre, width, delay)
% generate the combined gifs

r1 = median(icepat(:));
r2 = median(vicpaas(:));

vicpaas = vicpaas * r1/r2;

combinedData = [icepat vicpaas];
headerCombined = header;
headerCombined.sizeX = size(combinedData, 2);
headerCombined.sizeY = size(combinedData, 1);

analyze2gifWithWindowSetting(combinedData, headerCombined, giffile, sizeRatio, centre, width, delay);
