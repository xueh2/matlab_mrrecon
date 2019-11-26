
function [data, header, acq_time, physio_time, dicom_header] = readGTExportImageSeries_CMR_View(folderName, seriesNum)
% read in the gt create images (after CMR view) and get dicom headers
% Usage:
% [data, header, acq_time, physio_time, dicom_header] = readGTExportImageSeries_CMR_View(folderName, seriesNum)

if(nargin==1)
    seriesNum = folderName;
    folderName = '.';
end


[data, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries(folderName, seriesNum, 1, 0);
data = double(data);

if(nargout>2)
    acq_time = squeeze(acq_time);
end

if(nargout>3)
    physio_time = squeeze(physio_time);
end

RO = size(data,1);
E1 = size(data,2);
SLC = size(data,3);
E2 = size(data,4);
CON = size(data,5);
PHS = size(data,6);
REP = size(data,7);
SET = size(data,8);
AVE = size(data,9);

dicom_header = [];

for slc=1:SLC
    im = data(:,:,slc,1,1,1,1,1,1,1);
    gt_h = header(slc,1,1,1,1,1,1);
    im = squeeze(im);
    im = permute(im, [2 1]);
    h = ComputeDicomCoordFromGtOffline(gt_h.PatientPosition, gt_h.read_dir, gt_h.phase_dir, double(gt_h.FOV([2 1 3])), [size(im) 1], size(im, 1), size(im, 2));
    dicom_header = [dicom_header; h];
end

data = squeeze(data);
header = squeeze(header);
