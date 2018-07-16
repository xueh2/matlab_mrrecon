
function plot_data_gaussian_ellipsoid(data1, data2, legend_str, axis_limits)
% plot_data_gaussian_ellipsoid(data, legend_str, axis_limits)

num_d = numel(data1);

for n=1:num_d
    disp(['data ' num2str(n) ' mean/std : ' num2str(mean(data1{n})) '+/-' num2str(std(data1{n})) ' - ' num2str(mean(data2{n})) '+/-' num2str(std(data2{n}))]);
end

SD = 1.0;
npts = 500;

deco = {'r+', 'b+', 'ko', 'gs', 'm^', 'c*', 'ks'};

figure;
hold on

ind = 1;
for n=1:num_d
    plot(data1{n}, data2{n}, deco{ind}, 'MarkerSize', 10);
    ind = ind+1;
end

hold off
axis(axis_limits);
box on
grid on
legend(legend_str, 'FontSize', 14);

ind = 1;
for n=1:num_d
    h = gca;
    dd = [data1{n}, data2{n}];
    md = mean(dd, 1);
    C = cov(dd);
    hp = plot_gaussian_ellipsoid(md, C, SD, npts, h);
    c_str = deco{ind};
    set(hp,'color', c_str(1), 'LineWidth', 3); 
    ind = ind + 1;
end
