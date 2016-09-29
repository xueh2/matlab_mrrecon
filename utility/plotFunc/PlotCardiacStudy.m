function PlotCardiacStudy(dataBuf, indToPlot, rowCol)

dataAll = zeros([size(dataBuf{indToPlot(1),3}) numel(indToPlot)]);

for kk=1:numel(indToPlot)
    dataAll(:,:,:,kk) = dataBuf{indToPlot(kk), 3};
end

figure; imagescn(dataAll, [], rowCol, [], 3);