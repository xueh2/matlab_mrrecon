function [data, snr, gmap] = readSNRReconRes(home, dataSeriesNum, snrSeriesNum, gmapSeriesNum, plotFlag)

gmap = readGTPlusExportImageSeries(home, gmapSeriesNum);
gmap = squeeze(gmap);

if(plotFlag)
    figure; imagescn(gmap/100, [0 2], [], [], 4); colormap('jet')
end

snr = readGTPlusExportImageSeries(home, snrSeriesNum); snr = squeeze(snr);

data = readGTPlusExportImageSeries(home, dataSeriesNum); data = squeeze(data);