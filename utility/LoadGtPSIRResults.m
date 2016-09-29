
ori = readGTPlusExportImageSeries_Squeeze('.', 100);

ori_norm_MagIR = readGTPlusExportImageSeries_Squeeze('.', 101);
ori_MagPD = readGTPlusExportImageSeries_Squeeze('.', 102);
ori_PSIR = readGTPlusExportImageSeries_Squeeze('.', 103);

moco = readGTPlusExportImageSeries_Squeeze('.', 104);
moco_MagIR = readGTPlusExportImageSeries_Squeeze('.', 105);
moco_MagPD = readGTPlusExportImageSeries_Squeeze('.', 106);
moco_PSIR = readGTPlusExportImageSeries_Squeeze('.', 107);

moco_ave = readGTPlusExportImageSeries_Squeeze('.', 108);
moco_ave_norm_MagIR = readGTPlusExportImageSeries_Squeeze('.', 109);
moco_ave_MagPD = readGTPlusExportImageSeries_Squeeze('.', 110);
moco_ave_norm_PSIR = readGTPlusExportImageSeries_Squeeze('.', 111);

ori_MAGIR = readGTPlusExportImageSeries_Squeeze('.', 112);
moco_MAGIR = readGTPlusExportImageSeries_Squeeze('.', 113);
moco_ave_MAGIR = readGTPlusExportImageSeries_Squeeze('.', 114);