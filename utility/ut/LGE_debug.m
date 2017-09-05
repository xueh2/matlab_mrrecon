
cd E:\ut\LGE\VD13A_AERA_LGE_MOCO_BLACK_SPOT_CASE\meas_MID00059_FID94643_3minCV_gtPlus_PSIR_9sl_4ave_MOCO\moco_ave_res
cd E:\ut\LGE\meas_MID01986_FID19364_Sax_RadCV_gtPlus_PSIR_12sl_4ave_MOCO\moco_ave_res
cd E:\ut\LGE\MOCO_Evaluation\20150408_14h22m26s_73662\grappa_res

cd E:\gtuser\gt_windows_setup\ut\LGE\gtPSIR-2015-09-30\meas_MID00268_FID61977_PSIR_4ch_MOCO\moco_ave_res
cd \\hl-share\RawMRI\Lab-Kellman\Share\temp\DB_LGE_MOCO_AVE_OnTheFly_42363_13415235_13415240_75_20170728-151249\moco_ave_res

% early enhancement
cd E:\ut\LGE\MOCO_Evaluation\20150408_14h04m58s_73641\grappa_res
cd E:\ut\LGE\MOCO_Evaluation\20150408_14h02m51s_73640\grappa_res

% moco_ave
moco_ave_MAGIR = readGTPlusExportImageSeries_Squeeze('.', 109);
moco_ave_PSIR = readGTPlusExportImageSeries_Squeeze('.', 111);
moco_ave_PD = readGTPlusExportImageSeries_Squeeze('.', 110);
moco_ave_MAGIR_noSCC = readGTPlusExportImageSeries_Squeeze('.', 114);

figure; imagescn(cat(4, moco_ave_MAGIR, moco_ave_PSIR), [], [2 size(moco_ave_MAGIR, 3)], 15);

figure; imagescn(cat(4, moco_ave_MAGIR, moco_ave_PSIR, moco_ave_PD, moco_ave_MAGIR_noSCC), [], [4 9], 15);

s = size(moco_ave_PSIR);
moco_ave_PSIR = Matlab_gt_resize_2D_image(moco_ave_PSIR, 4*s(1), 4*s(2), 5);

figure; imagescn( moco_ave_PSIR(:,:, end:-1:1), [], [3 3]);

[I, xlim, ylim, clim, position]=getimagedata(gcf);
setimagescn(xlim,ylim,clim);

moco_ave = readGTPlusExportImageSeries_Squeeze('.', 108);

moco_ave = readGTPlusExportImageSeries_Squeeze('.', 108);
figure; imagescn( abs(moco_ave), [], [], 15, 4);

d = moco_ave(:,:,6, :);
d = squeeze(d);
figure; imagescn(abs(d));

PD = d(:,:,2);
imagescn(angle(PD))

a = PD./abs(PD);

IR = d(:,:,1);
PSIR = IR .* conj(a);
imagescn(real(PSIR))

% ori
ori = readGTPlusExportImageSeries_Squeeze('.', 100);
figure; imagescn(ori, [], [2 size(ori, 3)], 15, 5);

ori_MAGIR = readGTPlusExportImageSeries_Squeeze('.', 101);
ori_PSIR = readGTPlusExportImageSeries_Squeeze('.', 103);
ori_PD = readGTPlusExportImageSeries_Squeeze('.', 102);
ori_MAGIR_noSCC = readGTPlusExportImageSeries_Squeeze('.', 112);

figure; imagescn(cat(4, ori_MAGIR, ori_MAGIR_noSCC), [], [2 9], 15, 5);