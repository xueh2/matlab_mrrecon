function plotTwoAIFLut(params, T1, T2, T1_2, T2_2)
% plot two LUT for AIF
% plotTwoAIFLut(params, T1, T2, T1_2, T2_2)

params.T1_0 = T1/1e3;
params.T2_0 = T2/1e3;
[LUT1, Gd1] = perfusion_lut_flash_pd_flash_sr(params);
l1 = ['T1_0=' num2str(T1) 'ms; T2_0=' num2str(T2) 'ms']

params.T1_0 = T1_2/1e3;
params.T2_0 = T2_2/1e3;
[LUT2, Gd2] = perfusion_lut_flash_pd_flash_sr(params);
l2 = ['T1_0=' num2str(T1_2) 'ms; T2_0=' num2str(T2_2) 'ms']

figure
hold on
plot(Gd1, LUT1);
plot(Gd2, LUT2, 'r');
hold off
legend(l1, l2);
box on
grid on
title('AIF, Gd LUT');
xlabel('Gd, mmol/l');
ylabel('SR/PD')