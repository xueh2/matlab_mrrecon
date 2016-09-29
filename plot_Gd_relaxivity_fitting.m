function r = plot_Gd_relaxivity_fitting(Gd, tx, contrast, RxStr)
% plot Gd r1/r2 fitting
% plot_Gd_relaxivity_fitting(Gd, tx, contrast)

tx_used = tx(1:numel(Gd));
Rx_used = 1./tx_used;

% Rx_used = Rx_used - Rx_used(1);

P = polyfit(Gd(:), 1000*Rx_used(:),1);
Rx_fit = P(1)*Gd+P(2);

xpos = 1
ypos = (Rx_fit(end-1) + Rx_fit(end-2)) / 2
theString = sprintf('y = %.5f x + %.5f', P(1), P(2));

figure; 
hold on
plot(Gd(:), 1000*Rx_used(:), '+');
plot(Gd(:), Rx_fit(:),'r-.');
text(xpos, ypos, theString, 'FontSize', 12);
hold off
xlabel('Gd, mmol/l');
ylabel([RxStr '-' RxStr '_0, Hz']);
title(contrast)
box on

r = P(1);