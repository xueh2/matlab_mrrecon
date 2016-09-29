function [ori, ori_psir, ori_mag_PD, ori_mag_ir] = ReviewLGEPSIR(resDir)
% ReviewLGEPSIR(resDir)

ori = readGTPlusExportImageSeries_Squeeze(resDir, 100);
ori_psir = readGTPlusExportImageSeries_Squeeze(resDir, 103);
ori_mag_PD = readGTPlusExportImageSeries_Squeeze(resDir, 102);
ori_mag_ir = readGTPlusExportImageSeries_Squeeze(resDir, 112);

scalingFactor = 8;

figure('Name','ORI MAGIR','NumberTitle','off');
imagescn(ori_mag_ir, [], [], scalingFactor);

figure('Name','ORI PSIR','NumberTitle','off');
imagescn(ori_psir, [], [], scalingFactor);

figure('Name','ORI MAG PD','NumberTitle','off');
imagescn(ori_mag_PD, [], [], scalingFactor);

figure('Name','ORI','NumberTitle','off');
imagescn(ori, [], [], scalingFactor);
