function plot_x_y_with_linefit(x, y, x_label, y_label, title_str)
% plot_x_y_with_linefit(x, y, x_label, y_label)

figure
hold on
plot(x, y, '.');
lsline
P = polyfit(x, y, 1);
eqn = ['y = ' sprintf('%3.5fx + %3.5f', P(1), P(2))];
hold off
box on
grid on
xlabel(['x : ' x_label])
xlabel(['y : ' y_label])
axis equal
lx = get(gca, 'XLim')
ly = get(gca, 'YLim')
text(lx(1) + 0.05*(lx(2)-lx(1)), ly(2)*0.85,[ eqn ])
title(title_str)
