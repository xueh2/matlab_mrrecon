
function h = PlotContourForAHAModel(AHA_s1, h)
% h = PlotContourForAHAModel(AHA_s1, h)

num_sectors = numel(AHA_s1);

cmap = hsv(num_sectors);

if(ishandle(h))
    axes(h);
else
    figure;
    h = gca;
end

hold on
for ii=1:num_sectors
    C3 = AHA_s1{ii};
    plot(h, C3(:,1), C3(:,2), '-', 'Color',cmap(ii,:));
end 
hold off