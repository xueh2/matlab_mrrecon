
function [data, header, acq_time, physio_time, attribs, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries(folderName, seriesNum, withTime, numAsRep)
% read in the gtplus create images
% data = readGTPlusExportImageSeries(folderName, seriesNum);
% [data, header, acq_time, physio_time, endo_pt, epi_pt] = readGTPlusExportImageSeries(folderName, seriesNum, withTime, numAsRep)

if nargin < 3
    withTime = 0;
end

if nargin < 4
    numAsRep = 0;
end

[names, num] = findFILE(folderName, '*.hdr');

endo_pt = [];
epi_pt = [];
user_int = [];
attribs = [];

output_format = getenv('OutputFormat');

if(num==0 | strcmp(output_format, 'h5'))
    [data, header, attribs, acq_time, physio_time] = readGTPlusExportImageSeries_h5(folderName, seriesNum, withTime, numAsRep);
else
    header = [];
    [data, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_hdr(folderName, seriesNum, withTime, numAsRep);
end



