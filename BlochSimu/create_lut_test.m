
params.r1 = 4.9188;
params.r2 = 5.9133;
params.T1_0 = 1.2;
params.T2_0 = 0.045;

%% SR_SSFP & PD_FLASH LUT

params.TR = 0.00243;
params.TD = 0.0376;
params.PD_flip_angle = 5;
params.SR_flip_angle = 49.9999;
params.offresonance = 0;
params.Npe_full = 111;
params.PAT_accel_factor =  3;
params.steadystate_prep = 'linear';
params.N_runup = 3;
params.rf_phase_spoiler_increment = 112;

[LUT, Gd] = perfusion_lut_flash_pd_ssfp_sr(params);
LUT_gt = Matlab_gt_perfusion_bloch_simulation(single(Gd'), params.PD_flip_angle, params.SR_flip_angle, params.TD*1e3, params.TR*1e3, params.Npe_full, params.PAT_accel_factor, 'SSFP', 'Flash', params.T1_0*1e3, params.T2_0*1e3, params.r1, params.r2, [] );
figure; hold on; plot(Gd, LUT); plot(Gd, LUT_gt, 'r+'); hold off; grid on

cd E:\gtuser\mrprogs\gadgetron\toolboxes\gtplus\ut\result
a = analyze75read('lut_flash_pd_ssfp_perf.hdr');
figure; hold on; plot(Gd, LUT); plot(Gd, a, 'r.'); hold off; grid on

%% SR_FLASH & PD_FLASH LUT
params.TR = 0.00243;
params.TD = 0.0376;
params.PD_flip_angle = 5;
params.SR_flip_angle = 14;
params.offresonance = 0;
params.Npe_full = 111;
params.PAT_accel_factor = 3;
params.steadystate_prep = 'linear';
params.N_runup = 3;
params.rf_phase_spoiler_increment = 112;

[LUT, Gd] = perfusion_lut_flash_pd_flash_sr(params);
figure; plot(Gd, LUT); grid on

LUT_gt = Matlab_gt_perfusion_bloch_simulation(single(Gd'), params.PD_flip_angle, params.SR_flip_angle, params.TD*1e3, params.TR*1e3, params.Npe_full, params.PAT_accel_factor, 'Flash', 'Flash', params.T1_0*1e3, params.T2_0*1e3, params.r1, params.r2, [] );
figure; hold on; plot(Gd, LUT); plot(Gd, LUT_gt, 'r+'); hold off; grid on

cd E:\gtuser\mrprogs\gadgetron\toolboxes\gtplus\ut\result
a = analyze75read('lut_flash_pd_flash_perf.hdr');
figure; hold on; plot(Gd, LUT); plot(Gd, a, 'r.'); hold off; grid on

%% SR_FISP& PD_FISP LUT
params.TR = 0.00243;
params.TD = 0.0376;
params.PD_flip_angle = 5;
params.SR_flip_angle = 50;
params.offresonance = 0;
params.Npe_full = 111;
params.PAT_accel_factor = 3;
params.steadystate_prep = 'linear';
params.N_runup = 3;
params.rf_phase_spoiler_increment = 112;

[LUT, Gd] = perfusion_lut_ssfp_pd_ssfp_sr(params);
figure; plot(Gd, LUT); grid on

LUT_gt = Matlab_gt_perfusion_bloch_simulation(single(Gd'), params.PD_flip_angle, params.SR_flip_angle, params.TD*1e3, params.TR*1e3, params.Npe_full, params.PAT_accel_factor, 'SSFP', 'SSFP', params.T1_0*1e3, params.T2_0*1e3, params.r1, params.r2, [] );
figure; hold on; plot(Gd, LUT); plot(Gd, LUT_gt, 'r+'); hold off; grid on

cd E:\gtuser\mrprogs\gadgetron\toolboxes\gtplus\ut\result
a = analyze75read('lut_ssfp_pd_ssfp_perf.hdr');
figure; hold on; plot(Gd, LUT); plot(Gd, a, 'r.'); hold off; grid on

%% AIF LUT
% params.TR = 0.00112; % 1.12 ms (1-echo AIF protocol)
params.TR = 0.00212; % 2.12 ms (2 echo AIF protocol)
params.TD = 0.005;
params.PD_flip_angle = 8;
params.SR_flip_angle = 8;
params.offresonance = 0;
params.Npe_full = 36;
params.PAT_accel_factor = 2;
params.steadystate_prep = 'linear';
params.N_runup = 3;
params.rf_phase_spoiler_increment = 112;

[LUT, Gd] = perfusion_lut_flash_pd_flash_sr(params);
figure; plot(Gd, LUT); grid on

LUT_gt = Matlab_gt_perfusion_bloch_simulation(single(Gd'), params.PD_flip_angle, params.SR_flip_angle, params.TD*1e3, params.TR*1e3, params.Npe_full, params.PAT_accel_factor, 'Flash', 'Flash', params.T1_0*1e3, params.T2_0*1e3, params.r1, params.r2, [] );
figure; hold on; plot(Gd, LUT); plot(Gd, LUT_gt, 'r+'); hold off; grid on

cd E:\gtuser\mrprogs\gadgetron\toolboxes\gtplus\ut\result
a = analyze75read('AIF_lut_flash_pd_flash_perf.hdr');
figure; hold on; plot(Gd, LUT); plot(Gd, a, 'r.'); hold off; grid on

%% test the gt bloch simulation
cd D:\gtuser\gt_windows_setup\ut\DualBolus\20150506\stress\meas_MID00300_FID45790_BIG_STRESS__075_PK_FLASH_PERFUSION_GT_AIF2eTPAT2\DebutOutput

SRNorm = analyze75read('SRNorm_0.hdr');
size(SRNorm)
SRNorm = single(SRNorm);

% SR_FLASH & PD_FLASH LUT
params.TR = 0.0027;
params.TD = 0.005;
params.PD_flip_angle = 5;
params.SR_flip_angle = 14;
params.offresonance = 0;
params.Npe_full = 111;
params.PAT_accel_factor = 3;
params.steadystate_prep = 'linear';
params.N_runup = 5;
params.rf_phase_spoiler_increment = 112;

[LUT, Gd] = perfusion_lut_flash_pd_flash_sr(params);

LUT_gt = Matlab_gt_perfusion_bloch_simulation(single(Gd'), params.PD_flip_angle, params.SR_flip_angle, params.TD*1e3, params.TR*1e3, params.Npe_full, params.PAT_accel_factor, 'Flash', 'Flash', [] );
figure; hold on; plot(Gd, LUT); plot(Gd, LUT_gt, 'r+'); hold off; grid on

[LUT_gt, CASignal_gt] = Matlab_gt_perfusion_bloch_simulation(single(Gd'), params.PD_flip_angle, params.SR_flip_angle, params.TD*1e3, params.TR*1e3, params.Npe_full, params.PAT_accel_factor, 'Flash', 'Flash', SRNorm );
CASignal = interp1(LUT_gt, Gd, SRNorm, 'linear');

figure; imagescn(cat(4, SRNorm, CASignal, CASignal_gt), [], [], [], 3);

LUT_gt2 = analyze75read('Perf_T1_Correction_LUT.hdr');
CASignal2 = interp1(LUT_gt2, Gd, SRNorm, 'linear');
CASignal_gt2 = analyze75read('CASignal_LUTCorrection_0.hdr');
figure; imagescn(cat(4, SRNorm, CASignal2, CASignal_gt2), [], [], [], 3);

aif_echo0 = analyze75read('aif_cin_all_echo0_signal.hdr');
aif_echo1 = analyze75read('aif_cin_all_echo1_signal.hdr');
figure; hold on; plot(aif_echo0); plot(aif_echo1, 'r.'); hold off; grid on

aif_R2Star = analyze75read('aif_cin_all_R2Star.hdr');
aif_R2Star_SLEP = analyze75read('aif_cin_all_R2Star_SLEP.hdr');
figure; hold on; plot(aif_R2Star); plot(aif_R2Star_SLEP, 'r.'); hold off; grid on

aif_cin_all_echo0_signal_after_R2StarCorrection = analyze75read('aif_cin_all_echo0_signal_after_R2StarCorrection.hdr');
figure; 
hold on; 
plot(aif_echo0); 
plot(aif_echo1, 'r.'); 
plot(aif_cin_all_echo0_signal_after_R2StarCorrection, 'k'); 
hold off; 
grid on
legend('aif echo0', 'aif echo1', 'aif echo0 with T2* correction');

foot = 34;
peak = 52;

aif_LUT_gt = analyze75read('aif_cin_LUT_flash_pd_flash_sr.hdr');
aif_cin_all_echo0_signal_after_R2StarCorrection_divide_PD_echo0 = analyze75read('aif_cin_all_echo0_signal_after_R2StarCorrection_divide_PD_echo0.hdr');
aif_cin_all_echo0_LUTCorrection = analyze75read('aif_cin_all_echo0_LUTCorrection.hdr');
aif_CASignal = interp1(single(aif_LUT_gt), Gd, aif_cin_all_echo0_signal_after_R2StarCorrection_divide_PD_echo0, 'linear');
figure; hold on; plot(aif_CASignal); plot(aif_cin_all_echo0_LUTCorrection, 'r.'); hold off; grid on

aif_cin = analyze75read('aif_cin.hdr');
figure; plot(aif_cin);

perf_mask_0 = analyze75read('perf_mask_0.hdr');
perf_moco_0 = analyze75read('perf_moco_0.hdr');
perf_moco_upsampled_0 = analyze75read('perf_moco_upsampled_0.hdr');
figure; imagescn(cat(4, CASignal_gt2, perf_moco_0), [], [], [], 3);
figure; imagescn(perf_moco_upsampled_0, [], [], [], 3);

perf = perf_moco_upsampled_0(:,:,foot+1:end);
figure; imagescn(perf, [], [], [], 3);

x = 107;
y = 84;

perf_s = perf(x, y, :);
perf_s = squeeze(perf_s);
figure; hold on; plot(aif_cin); plot(perf_s, 'r.'); hold off; grid on

cin = single(aif_cin);
y = single(perf_s);
orderBSpline = 4;
numOfInternalControlPoints_L1 = 4;
lambda = 0.001;
[r, yr] = PerformDeconvolution_L1_Fista_BSpline(cin, y, orderBSpline, numOfInternalControlPoints_L1, lambda);

max_iter_L1BSpline = 100;
lambda_L1BSpline = 0.001;
obj_thres_L1BSpline = 1e-6;
grad_thres_L1BSpline = 1e-6;
print_iter_L1BSpline = 0;
num_of_wavLevels_L1BSpline = 1;
with_approx_coeff_L1BSpline = 1;
max_iter_Fermi = 100;
[r_gt, yr_gt, r_Fermi, yr_Fermi] = Matlab_gt_deconvolution(cin, y, 'L1BSpline', orderBSpline, numOfInternalControlPoints_L1, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi);

figure; hold on; plot(r); plot(r_gt, 'r.'); plot(r_Fermi, 'g.'); hold off; grid on

figure; hold on; plot(y); plot(yr, 'r.');  plot(yr_gt, 'k.'); plot(yr_Fermi, 'g.'); hold off; grid on

deltaT = 323.673;
data_length_FPWH_ratio = 3;
flowmaps = Matlab_gt_perfusion_flow_mapping(cin, single(perf_moco_upsampled_0), single(perf_mask_0), foot, peak, deltaT, data_length_FPWH_ratio, orderBSpline, numOfInternalControlPoints_L1, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi, []);
figure; imagescn(flowmaps);

perf_moco_0 = analyze75read('Input_perf_computeFlowMap_0.hdr');
cin = analyze75read('Input_cin_computeFlowMap.hdr');

deltaT = 323.673;
data_length_FPWH_ratio = 3;
flowmaps = Matlab_gt_perfusion_flow_mapping(cin, single(perf_moco_0), single(perf_mask_0), foot, peak, deltaT, data_length_FPWH_ratio, orderBSpline, numOfInternalControlPoints_L1, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi, []);
figure; imagescn(flowmaps);

cd D:\gtuser\mrprogs\install\DebugOutput
map0 = analyze75read('flow_maps_0.hdr');
map1 = analyze75read('flow_maps_1.hdr');
map2 = analyze75read('flow_maps_2.hdr');
figure; imagescn(cat(3, map0, map1, map2));

ms = map0;
ms = imresize(ms, 2);
mr = a(:,:,1,4)/1000;
figure; imagescn(cat(3, ms, mr)); PerfColorMap
