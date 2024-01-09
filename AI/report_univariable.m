function [min_v, max_v, median_v, v_0p05, v_0p95, mean_v, std_v] = report_univariable(v, desp)

ind = find(v>0 & ~isnan(v));

mean_v = mean(v(ind));
std_v = std(v(ind));
min_v = min(v(ind));
max_v = max(v(ind));
median_v = median(v(ind));
r = quantile(v(ind),[0.05 0.95]);
v_0p05 = r(1);
v_0p95 = r(2);

disp([desp ' - mean +/- std = ' num2str(mean_v) '+/-' num2str(std_v) ' - min = ' num2str(min(v(ind))) ' - max = ' num2str(max(v(ind))) ' - median = ' num2str(median(v(ind))) ' - 5% - 95% = ' num2str(prctile(v(ind), 5)) ' - ' num2str(prctile(v(ind), 95))]);
