
function [data1, data2] = compareGTPlusExportData(baseFileName1, baseFileName2)

data1 = readGTPlusExportData(baseFileName1);
data2 = readGTPlusExportData(baseFileName2);

norm(data1(:))
norm(data2(:))
norm(data1(:)-data2(:))