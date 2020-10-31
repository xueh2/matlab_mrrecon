
function [data, header, acq_time, physio_time, dicom_header] = readGTExportImageSeries_without_CMR_View(folderName, seriesNum)
% read in the gt create images (before CMR view) and get dicom headers in cmr view
% Usage:
% [data, header, acq_time, physio_time, dicom_header] = readGTExportImageSeries_without_CMR_View(folderName, seriesNum)

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
    im = data(:,:,slc,1,CON,PHS,REP,SET,AVE,1);
    im = squeeze(im);
    if(numel(size(header))==3)
        gt_h = header(slc,1,end);
    end
    if(numel(size(header))==4)
        gt_h = header(slc,1,end,end);
    end
    if(numel(size(header))==5)
        gt_h = header(slc,1,end,end,end);
    end
    if(numel(size(header))==6)
        gt_h = header(slc,1,end,end,end,end);
    end
    if(numel(size(header))==7)
        gt_h = header(slc,1,end,end,end,end,end);
    end
    
    im = flipdim(flipdim(im, 2), 1);

    ImageOrientationPatient = [gt_h.read_dir gt_h.phase_dir gt_h.slice_dir];
    [rotate90, flip] = NormOrientation_rev3 (ImageOrientationPatient)
    [I, gt_h_cmr] = apply_normorientation(im, rotate90, flip, gt_h);
    
    h = ComputeDicomCoordFromGtOffline(gt_h.PatientPosition, gt_h.read_dir, gt_h.phase_dir, double(gt_h.FOV), [size(im) 1], size(im, 1), size(im, 2));
    dicom_header = [dicom_header; h];
end

data = squeeze(data);
header = squeeze(header);
