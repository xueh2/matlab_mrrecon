

clear all
close all

%% compare the unit test results with ground truth

GTHome = getenv('GADGETRON_HOME')
addpath(fullfile(GTHome, 'bin'));
addpath(fullfile(GTHome, 'lib'));

UTDir = getenv('GTPLUS_UT_DIR')
addpath(UTDir)

cd(fullfile(GTHome, 'DebugOutput'));

data = readGTPlusExportData('fullkspace_');

data = readGTPlusExportData('fullkspace_After');

data = readGTPlusExportData('kspace_before_PF_Filter');
data = readGTPlusExportData('kspace_after_PF_Filter');
data = readGTPlusExportData('kspace_after_POCS');
data = readGTPlusExportData('ComplexIm_afterRefFill');
data = readGTPlusExportData('coilMap_fullres');
data = readGTPlusExportData('aveComplexIm');

pp = reshape(data, [256 256 8 4*31]);
plotKSpaceArray(pp);
plotKSpaceSamplingPattern(pp(:,:,2,:));

plotComplexImageArray(pp);

pp2 = pp;
pp2(:, 141:end, :, :) = 0;
plotKSpaceArray(pp2);

plotKSpaceSamplingPattern(pp2(:,:,2,:));

data = readGTPlusExportData('incomingKSpace');
ref = readGTPlusExportData('incomingRef');

ref = squeeze(ref);
size(ref)

ref2 = ref(:,117:140,:,:,:);

plotKSpaceSamplingPattern(ref2(:,:,2,:, 1));

ImRef = SoS_TemporalArray(ref(:,:,:,:,1));
imagescn(ImRef);

data = squeeze(data);
size(data)

plotKSpaceSamplingPattern(data(:,:,2,:, 1));

Im = SoS_TemporalArray(data(:,:,:,:,1));

imagescn(Im);

coilMap = readGTPlusExportData('coilMap_');
imagescn(squeeze(abs(coilMap)), [], [], [], 4);

src = readGTPlusExportData('ref_src_');
dst = readGTPlusExportData('ref_dst_');

refCoilMap = readGTPlusExportData('refCoilMap_filtered');

src = ref2(:,:,:,12, 1);
dst = src;

headerSrc = CreateFtkHeaderInfo(src, [1 1 1]);
headerDst = CreateFtkHeaderInfo(dst, [1 1 1]);
[ker, kerIm, unmixing, gFactor] = Matlab_PerformVICPAAGKernelCalibrationSrcDst2D(single(src), headerSrc, single(dst), headerDst, 5, [-4 0 4 8], 1, 0, [256 256], 0.0001, 4, 1);

csm = CoilSensitivity_Souheil_Walsh(ifft2c(refCoilMap), [7 7], 0);
imagescn(abs(csm));

aliasedIm = ifft2c(data(:,:,:,12,:));
aliasedIm = squeeze(aliasedIm);
size(aliasedIm)

im = sum(aliasedIm(:,:,:,1).*unmixing, 3);
im = squeeze(im);
imagescn(abs(im))

cd D:\gtuser\gt_windows_setup\ut\HASTE\raw30488
cd D:\gtuser\gt_windows_setup\ut\cine\meas_MID00832_FID37686_PK_realtime_tfi_9slice_tpat5_724
cd D:\gtuser\gt_windows_setup\ut\cine\meas_MID00832_FID37686_PK_realtime_tfi_9slice_tpat5_724\grappa_ref

cd D:\gtuser\gt_windows_setup\ut\cine\meas_MID02928_FID20308_PK_CV_gtPlus_RealTimeCine_NonLinear

data = readGTPlusExportImages('.', 'GT_Image');
size(data)

data = squeeze(data);
header = CreateGtImageHeader(data(:,:,:,1), [0 0 0], [1.09 1.09 3]);

for i=1:4
    v1 = data(:,:,:,i);
    v2 = v1;
    v2(:,:,1:2:31) = v1(:,:,1:16);
    v2(:,:,2:2:31) = v1(:,:,1:15);
    Matlab_gt_write_analyze(single(v2), header, ['volume' num2str(i)]);
end

moco = readGTPlusExportImages('.', 'GT_Image_MOCO');
moco = squeeze(moco);
size(moco)

s = std(data, 0, 4);

