
function [perf_corrected, LUT] = Perfusion_PrepareSignal_Perf_T1Correction(perf_SROverPD, TE, FA_SR, FA_PD, TR, TS, TD, E1_full, accel_factor, seq_type, seq_type_PD, T1_0, T2_0, r1, r2, LUTCorrection)
% perform the perf T1 correction
% all times are in ms
% TE: echo time
% FA_SR, FA_PD : flip angle for AIF SR and PD images
% TR, TS, TD : aif readout TR, saturation time (from the end of SR pulse to the center kspace line) and TD (from end of SR pulse to the binning of first readout)
% E1_full, acce_factor
% seq_type : SSFP or Flash
% Gd_method : TwoPoint or LUT

%% compute mean PD

RO = size(perf_SROverPD, 1);
E1 = size(perf_SROverPD, 2);

REP = size(perf_SROverPD, 3);

% T1 correction
if(LUTCorrection)
    Gd = [0:0.01:20]';     
    [LUT, SR_Gd] = Matlab_gt_perfusion_bloch_simulation(Gd, FA_PD, FA_SR, TD, TR, E1_full, accel_factor, seq_type, seq_type_PD, T1_0, T2_0, r1, r2, single(perf_SROverPD));
else
    ratio = sin(FA_PD*pi/180)/sin(FA_SR*pi/180);
    SR_Gd = perf_SROverPD * ratio;
    SR_Gd = log(1-SR_Gd) * (-1/TS);
      
end

perf_corrected = SR_Gd;
