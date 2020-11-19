function [min_v, max_v, median_v, v_0p25, v_0p75] = report_univariable(v, desp)

ind = find(v>0);

mean_v = mean(v(ind));
std_v = std(v(ind));
min_v = min(v(ind));
max_v = max(v(ind));
median_v = median(v(ind));
r = quantile(v(ind),[0.25 0.75]);
v_0p25 = r(1);
v_0p75 = r(2);

disp([desp ' - mean +/- std = ' num2str(mean_v) '+/-' num2str(std_v) ' - min = ' num2str(min(v(ind))) ' - max = ' num2str(max(v(ind))) ' - median = ' num2str(median(v(ind))) ' - 25% - 75% = ' num2str(quantile(v(ind),[0.25 0.75]))]);
