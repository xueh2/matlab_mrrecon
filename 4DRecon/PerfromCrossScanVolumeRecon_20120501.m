
% perform the cross scan volume recon fusion
pixelSpacing = [1.82 1.82 1.82];

N = numel(dataDirs)
volumeLinearAve = cell(N, 1);
volumeMoCoAve = cell(N, 1);
volumeImagePositionPatient = cell(N, 1);
volumeImageOrientationPatient = cell(N, 1);
volumeMoCoAveBFFD = cell(N, 1);
headerMoCoAveBFFD = cell(N, 1);

for n=1:N
    dataDir = dataDirs{n}
    cd(dataDir)

    load(fullfile(dataDir, option.mocoDir, 'AfterBreathingRetrogatingResults.mat'));
    load(fullfile(dataDir, option.mocoDir, 'MoCOAveBFFD.mat'));

    volumeLinearAve{n} = binAveAll_BreathingGating_Linear;
    volumeMoCoAve{n} = binAveAll_BreathingGating_MocoAve;
    volumeImagePositionPatient{n} = binAveAll_BreathingGating_ImagePositionPatient;
    volumeImageOrientationPatient{n} = binAveAll_BreathingGating_ImageOrientationPatient;
    volumeMoCoAveBFFD{n} = volumeBFFD;
    headerMoCoAveBFFD{n} = headerBFFD;
end

NPhs = size(volumeMoCoAveBFFD{1}, 4);

imageSize = size(binAveAll_BreathingGating_Linear);

filename = fullfile(dataDirs{1}, option.mocoDir, 'FFDMask.mat');
if ( ~isFileExist(filename) )
    [volumeMasked, headerMasked, mask] = PerformMaskFFDVolumeWithAcquiredRegion(volumeMoCoAveBFFD{1}(:,:,:,1), headerMoCoAveBFFD{1}, imageSize, volumeImagePositionPatient{1}, volumeImageOrientationPatient{1}, 'linear');
    save(filename, 'mask');
else
    load(filename);
end

volume1_mask = mask;
headerVolume1Mask = headerMoCoAveBFFD{1};
headerVolume1Mask.sizeT = 1;
Matlab_SaveAnalyze(single(mask), headerVolume1Mask, fullfile(dataDirs{1}, option.mocoDir, 'volume1_mask.hdr'));

ind = find(mask==0);
volumeLinear = volumeMoCoAveBFFD{1};
for phs = 1:NPhs
    aVolume = volumeLinear(:,:,:,phs);
    aVolume(ind(:)) = 0;
    volumeLinear(:,:,:,phs) = aVolume;
end
volumeMoCoAveBFFD{1} = volumeLinear;

Matlab_SaveAnalyze(single(volumeMoCoAveBFFD{1}), headerMoCoAveBFFD{1}, fullfile(dataDirs{1}, option.mocoDir, 'volume1.hdr'));

% [headers, imagePositionPatient2D, imageOrientationPatient2D] = CreateHeader2DFrameOf3DVolume(volumeMasked, headerMasked);
% Render4DVolume(volumeMasked, pixelSpacing, imagePositionPatient2D, imageOrientationPatient2D, -1, 1, 1, option.centre, option.width, 1/20, 0, []);
% headerMasked2 = headerMasked;
% headerMasked2.sizeT = 1;
% Matlab_SaveAnalyze(single(volumeMasked(:,:,:,1)), headerMasked2, 'volumeMasked1.hdr');

% --------------------

filename = fullfile(dataDirs{2}, option.mocoDir, 'FFDMask.mat');
if ( ~isFileExist(filename) )
    [volumeMasked, headerMasked, mask] = PerformMaskFFDVolumeWithAcquiredRegion(volumeMoCoAveBFFD{2}(:,:,:,1), headerMoCoAveBFFD{2}, imageSize, volumeImagePositionPatient{2}, volumeImageOrientationPatient{2}, 'rotate');
    save(filename, 'mask');
else
    load(filename);
end

volume2_mask = mask;
headerVolume2Mask = headerMoCoAveBFFD{2};
headerVolume2Mask.sizeT = 1;
Matlab_SaveAnalyze(single(mask), headerVolume2Mask, fullfile(dataDirs{2}, option.mocoDir, 'volume2_mask.hdr'));

ind = find(mask(:,:,:,1)==0);
volumeRotate = volumeMoCoAveBFFD{2};
for phs = 1:NPhs
    aVolume = volumeRotate(:,:,:,phs);
    aVolume(ind(:)) = 0;
    volumeRotate(:,:,:,phs) = aVolume;
end
volumeMoCoAveBFFD{2} = volumeRotate;

Matlab_SaveAnalyze(single(volumeMoCoAveBFFD{2}), headerMoCoAveBFFD{2}, fullfile(dataDirs{2}, option.mocoDir, 'volume2.hdr'));

% ---------------------
if ( ~isempty(xRange) )
    leftup = [xRange(1) yRange(1) zRange(1)];
    rightdown = [xRange(2) yRange(2) zRange(2)];

    [volume1_maskROI, headervolume1_maskROI]=getROI(volume1_mask, headerVolume1Mask, leftup, rightdown);
    Matlab_SaveAnalyze(int16(volume1_maskROI), headervolume1_maskROI, fullfile(dataDirs{1}, option.mocoDir, 'volume1_maskROI.hdr'));

else
    volume1_maskROI = volume1_mask;
    headervolume1_maskROI = headerVolume1Mask;
    Matlab_SaveAnalyze(int16(volume1_maskROI), headervolume1_maskROI, fullfile(dataDirs{1}, option.mocoDir, 'volume1_maskROI.hdr'));
end

% headerMasked2 = headerMasked;
% headerMasked2.sizeT = 1;
% Matlab_SaveAnalyze(single(volumeMasked(:,:,:,1)), headerMasked2, 'volumeMasked2.hdr');
% [headers, imagePositionPatient2D, imageOrientationPatient2D] = CreateHeader2DFrameOf3DVolume(volumeMasked, headerMasked);
% Render4DVolume(volumeMasked, pixelSpacing, imagePositionPatient2D, imageOrientationPatient2D, -1, 1, 1, option.centre, option.width, 1/20, 0, []);

% plot two scans

combinedFolder = fullfile(dataDirs{1}, option.mocoDir)

for phs=1:NPhs

cd(combinedFolder)

volume1 = volumeMoCoAveBFFD{1};
volume1(find(volume1<0)) = 0;
volume1 = volume1(:,:,:,phs);
size(volume1)
header1 = headerMoCoAveBFFD{1};
header1.sizeT = 1;
volume1 = volume1 * 10;

volume2 = volumeMoCoAveBFFD{2};
volume2(find(volume2<0)) = 0;
volume2 = volume2(:,:,:,phs);
size(volume2)
header2 = headerMoCoAveBFFD{2};
header2.sizeT = 1;
volume2 = volume2 * 10;

renderStepSize = 2;

[headers, imagePositionPatient2D, imageOrientationPatient2D] = CreateHeader2DFrameOf3DVolume(volume1, header1);
[headers2, imagePositionPatient2D2, imageOrientationPatient2D2] = CreateHeader2DFrameOf3DVolume(volume2, header2);

% Render4DVolume(volume1, pixelSpacing, imagePositionPatient2D, imageOrientationPatient2D, -1, phs, renderStepSize, option.centre, option.width, 1/20, 0, []);
% Render4DVolume(volume2, pixelSpacing, imagePositionPatient2D2, imageOrientationPatient2D2, -1, phs, renderStepSize, 200, 400, 1/20, 0, []);

% h = plotMrFtkImageIn3D(volume1(:,:,102), headers{102}, -1, 200, 400);
% h = plotMrFtkImageIn3D(volume2(:,:,102), headers2{102}, h, 200, 400);

%% resample one volume to the other's frame
[volumeResampled, headerResampled] = resampleVolume(volume2, header2, volume1, header1);

volumeResampled(find(volumeResampled<0)) = 0;
% Matlab_SaveAnalyze(int16(volume1), header1, 'volume1.hdr');
Matlab_SaveAnalyze(int16(volumeResampled), headerResampled, fullfile(combinedFolder, 'volumeResampled.hdr'));

%% cut the heart
% xRange = [40 150];
% yRange = [45 145];
% zRange = [55 175];

if ( ~isempty(xRange) )
    leftup = [xRange(1) yRange(1) zRange(1)];
    rightdown = [xRange(2) yRange(2) zRange(2)];

    [volumeResampledROI, headerResampledROI]=getROI(volumeResampled, headerResampled, leftup, rightdown);
    Matlab_SaveAnalyze(int16(volumeResampledROI), headerResampledROI, 'volumeResampledROI.hdr');

    [volume1ROI, header1ROI]=getROI(volume1, header1, leftup, rightdown);
    Matlab_SaveAnalyze(int16(volume1ROI), header1ROI, 'volume1ROI.hdr');
else
    volumeResampledROI = volumeResampled;
    headerResampledROI = headerResampled;
    Matlab_SaveAnalyze(int16(volumeResampledROI), headerResampledROI, 'volumeResampledROI.hdr');
    
    volume1ROI = volume1;
    header1ROI = header1;
    Matlab_SaveAnalyze(int16(volume1ROI), header1ROI, 'volume1ROI.hdr');
end

r = mean(volume1ROI(:))/mean(volumeResampledROI(:))

% volume1ROIMask = zeros(size(volume1ROI));
% volume1ROIMask(find(volume1ROI>0)) = 1;
% Matlab_SaveAnalyze(int16(volume1ROIMask), header1ROI, 'volume1ROIMask.hdr');
    
% get intensity scaling ratio
% xRange = [60 130];
% yRange = [65 125];
% zRange = [75 195];
% 
% leftup = [xRange(1) yRange(1) zRange(1)];
% rightdown = [xRange(2) yRange(2) zRange(2)];
% 
% [volumeResampledROI2, headerResampledROI2]=getROI(volumeResampled, headerResampled, leftup, rightdown);
% [volume1ROI2, header1ROI2]=getROI(volume1, header1, leftup, rightdown);
% 
% r = median(volume1ROI2(:))/median(volumeResampledROI2(:))

%% iterate the registration
command = ['D:\cc_views\hxue_pccs000608ws_snapshot_view2\SCR_MR_CV\MrFtk\prod\bin\vc10\Release\ftkRigidReg3D '...
    fullfile(combinedFolder, 'volumeResampledROI.hdr') ' ' fullfile(combinedFolder, 'volume1ROI.hdr') ...
    ' -parin D:\cc_views\hxue_pccs000608ws_snapshot_view2\SCR_MR_CV\MrFtk\ut\algorithm\registration\parameter\parameters3D.rigid -transform '...
    fullfile(combinedFolder, 'volume1ROIReg.hdr') ' -dofout ' fullfile(combinedFolder, 'volume1ROI_rigid.dof')]
dos(command, '-echo');

% transform the volume1ROIMask
command = ['D:\cc_views\hxue_pccs000608ws_snapshot_view2\SCR_MR_CV\MrFtk\prod\bin\vc10\Release\ftkTransformation3D '...
    fullfile(combinedFolder, 'volume1_maskROI.hdr') ' ' fullfile(combinedFolder, 'volume1ROIMask_Rigid.hdr') ' ' fullfile(combinedFolder, 'volume1ROI_rigid.dof') ...
    ' -Interpolator NN -Tp 0 ']
dos(command, '-echo');
            
command = ['D:\cc_views\hxue_pccs000608ws_snapshot_view2\SCR_MR_CV\MrFtk\prod\bin\vc10\Release\ftkAngioMOCO   1 1  '...
    fullfile(combinedFolder, 'volumeResampledROI.hdr') ' '...
    fullfile(combinedFolder, 'volume1ROIReg.hdr') ...
    '  -frameOfReferenceSame true  -regStrategy  Direct -rigid false  -stepDivRigid  5  -subsamplings 6 1  1  2  2  2  2 -rigidInitialFlag false -nonRigid true -invFlag true  -iters  6 100 100 100 100 100 100  -sigma  6  -neighbor  1  -stepDiv  3  -algo  GLCC  -numOfRep  1 -adaptiveSigmaForNonRigid true   -minimalSigma  24 -intensityClamp false -volumePreserving false']
dos(command, '-echo');

[volume1ROIMask_Rigid, header3D] = Matlab_LoadAnalyze(fullfile(combinedFolder, 'volume1ROIMask_Rigid.hdr'));
[dx3D, header3D] = Matlab_LoadAnalyze(fullfile(combinedFolder, 'volume1ROIReg_dx.hdr.hdr'));
[dy3D, header3D] = Matlab_LoadAnalyze(fullfile(combinedFolder, 'volume1ROIReg_dy.hdr.hdr'));
[dz3D, header3D] = Matlab_LoadAnalyze(fullfile(combinedFolder, 'volume1ROIReg_dz.hdr.hdr'));

volume1ROIMask_NonRigid = Matlab_PerformWarpingImage3D(volume1ROIMask_Rigid, header3D, single(dx3D), single(dy3D), single(dz3D), 'NRR', 5, 1);

ind = find(volume1ROIMask_NonRigid(:)<1);
volume1ROIMask_NonRigid(ind(:)) = 0;
Matlab_SaveAnalyze(int16(volume1ROIMask_NonRigid), header3D, fullfile(combinedFolder, ['volume1ROIMask_NonRigid.hdr']));

indInValidPixel_volume1ROI = find(volume1ROIMask_NonRigid(:)<=0);

[volume1ROIReg, header] = Matlab_LoadAnalyze('volume1ROIReg_MOCO.hdr.hdr');

indROI = find(volumeResampledROI(:)<=0);
indROIReg = find(volume1ROIReg(:)<=0);

volume1ROIReg = volume1ROIMask_NonRigid.*volume1ROIReg;

volumeCombinedROI = (r*volumeResampledROI + volume1ROIReg)/2;
Matlab_SaveAnalyze(int16(volumeCombinedROI), headerResampledROI, fullfile(combinedFolder, ['volumeCombinedROI2_' num2str(phs) '.hdr']));

%% FFD based combination
Matlab_SaveAnalyze(single(volume1ROIReg), header, 'volume1ROIReg_Masked.hdr');

s = size(volumeResampledROI);
binData = zeros(s(1), s(2), 2*s(3), 1);
binData(:,:,1:s(3)) = r*volumeResampledROI;
binData(:,:,s(3)+1:end) = volume1ROIReg;

pixelSpacing = [headerResampledROI.spacingX headerResampledROI.spacingY headerResampledROI.spacingZ];

[headers, imagePositionPatient2D, imageOrientationPatient2D] = CreateHeader2DFrameOf3DVolume(volumeResampledROI, headerResampledROI);
imagePositionPatient2D = [imagePositionPatient2D; imagePositionPatient2D];
imageOrientationPatient2D = [imageOrientationPatient2D; imageOrientationPatient2D];
dstPixelSpacing = pixelSpacing;
numOfSubdivision = 9;
tmpFolder = 'd:/temp'

backgroundThres = 0;

[volumeBFFD, headerBFFD] = PerformScatterInterpolationFFD(binData, pixelSpacing, imagePositionPatient2D, imageOrientationPatient2D, dstPixelSpacing, numOfSubdivision, tmpFolder, backgroundThres);

Matlab_SaveAnalyze(int16(volumeBFFD), headerBFFD, fullfile(combinedFolder, ['volumeCombinedROI3_' num2str(phs) '.hdr']));

end

% load the volume combined results
filename = fullfile(combinedFolder, ['volumeCombinedROI2_' num2str(1) '.hdr']);
[volume, header] = Matlab_LoadAnalyze(filename);
volumeCombinedROIAll = zeros([size(volume) NPhs]);

for phs=1:NPhs
    filename = fullfile(combinedFolder, ['volumeCombinedROI2_' num2str(phs) '.hdr']);
    [volume, header] = Matlab_LoadAnalyze(filename);
    volumeCombinedROIAll(:,:,:,phs) = volume;
end

lowR = 0;
highR = 4096;
mag = normalizeImage2Range(volumeCombinedROIAll, lowR, highR);
header.sizeT = NPhs;
Matlab_SaveAnalyze(single(volumeCombinedROIAll), header, fullfile(combinedFolder, 'volumeCombinedROIAll_ADD.hdr'));

Matlab_SaveAnalyze(single(mag), header, fullfile(combinedFolder, 'magWindowed.hdr'));
Save4DVolumeFromToHdrSeq(mag, header, 'D:\data\4DRecon_20111206\newVolunteers\Volunteer3\result\Series14_20', 'ADD');

% FFD
filename = fullfile(combinedFolder, ['volumeCombinedROI3_' num2str(1) '.hdr']);
[volume, header] = Matlab_LoadAnalyze(filename);
volumeCombinedROIAll = zeros([size(volume) NPhs]);

for phs=1:NPhs
    filename = fullfile(combinedFolder, ['volumeCombinedROI3_' num2str(phs) '.hdr']);
    [volume, header] = Matlab_LoadAnalyze(filename);
    volumeCombinedROIAll(:,:,:,phs) = volume;
end
header.sizeT = NPhs;
Matlab_SaveAnalyze(single(volumeCombinedROIAll), header, fullfile(combinedFolder, 'volumeCombinedROIAll_FFD.hdr'));

lowR = 0;
highR = 4096;
mag = normalizeImage2Range(volumeCombinedROIAll, lowR, highR);
Save4DVolumeFromToHdrSeq(mag, header, 'D:\data\4DRecon_20111206\newVolunteers\Volunteer3\result\Series14_20', 'FFD');

[volume, header] = Matlab_LoadAnalyze('D:\data\4DRecon_20111206\20120223\selectedSeries\series30\moco_newRecon_lessReg_MTD_0p75\volumeCombinedROIAll_MultiFFD.hdr');
lowR = 0;
highR = 4096;
mag = normalizeImage2Range(volume, lowR, highR);
Save4DVolumeFromToHdrSeq(mag, header, 'D:\data\4DRecon_20111206\presentations\20120712\results\OldVolunteer\Series30_32', 'FFD');

cd D:\data\4DRecon_20111206\newVolunteers\Volunteer4\Series20\moco_newRecon_lessReg_MTD_0p75
load MoCOAveBFFD.mat

lowR = 0;
highR = 4096;
mag = normalizeImage2Range(volumeBFFD, lowR, highR);
Save4DVolumeFromToHdrSeq(mag, headerBFFD, 'D:\data\4DRecon_20111206\presentations\20120712\results\Volunteer4\Series20', 'FFD');
Matlab_SaveAnalyze(single(volumeBFFD), headerBFFD, 'D:\data\4DRecon_20111206\presentations\20120712\results\Volunteer4\volume_FFD.hdr');

