
function [data, header, acq_time, physio_time] = readGTPlusExportImageSeries(folderName, seriesNum, withTime, numAsRep)
% read in the gtplus create images
% data = readGTPlusExportImageSeries(folderName, seriesNum);
% [data, header, acq_time, physio_time] = readGTPlusExportImageSeries(folderName, seriesNum, withTime, numAsRep);

if nargin < 3
    withTime = 0;
end

if nargin < 4
    numAsRep = 0;
end

[names, num] = findFILE(folderName, '*.hdr');

if(num==0)
    [data, header, acq_time, physio_time] = readGTPlusExportImageSeries_h5(folderName, seriesNum, withTime, numAsRep);
else
    header = [];
    [data, acq_time, physio_time] = readGTPlusExportImageSeries_hdr(folderName, seriesNum, withTime, numAsRep);
end



