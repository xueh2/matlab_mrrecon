function [mag, phs, acq_time, physio_time] = loadBinningFlowImages(imgDir, imgRole, phsRole)
% [mag, phs, acq_time, physio_time] = loadBinningFlowImages(imgDir)
% load binning flow results

if nargin < 2
    imgRole = 'GT_ImageRetro';
    phsRole = 'GT_Phase';   
end

[mag, acq_time, physio_time] = readGTPlusExportImages(imgDir, imgRole);
mag = squeeze(mag); size(mag)
[phs, acq_time, physio_time] = readGTPlusExportImages(imgDir, phsRole);
phs = squeeze(phs); size(phs)
acq_time = squeeze(acq_time);
physio_time = squeeze(physio_time);

