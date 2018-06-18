
function [data, header, acq_time, physio_time, endo_pt, epi_pt] = readGTPlusExportImageSeries_Squeeze(folderName, seriesNum, withTime, numAsRep)
% read in the gtplus create images
% Usage:
% data = readGTPlusExportImageSeries_Squeeze(folderName, seriesNum)
% [data, header, acq_time, physio_time, endo_pt, epi_pt] = readGTPlusExportImageSeries_Squeeze(folderName, seriesNum, withTime, numAsRep)

if(nargin==1)
    seriesNum = folderName;
    folderName = '.';
end

if(nargin<3)
    withTime = 0;
end

if(nargin<4)
    numAsRep = 0;
end

[data, header, acq_time, physio_time, endo_pt, epi_pt] = readGTPlusExportImageSeries(folderName, seriesNum, withTime, numAsRep);
size(data)
data = squeeze(data);
data = double(data);
size(data)

if(nargout>2)
    acq_time = squeeze(acq_time);
end

if(nargout>3)
    physio_time = squeeze(physio_time);
end

