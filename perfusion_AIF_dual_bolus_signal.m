function perfusion_AIF_dual_bolus_signal(AIF, mask, AIF_Full, mask_Full, params)
% get AIF signal for mini/full bolus

TE0 = params.TE0;
TE1 = params.TE1;
FA_PD = params.FA_PD;
FA_SR = params.FA_SR;
TS = params.TS;
TD = params.TD;
TR = params.TR;
numPD = params.numPD;
percentageAIF = params.percentageAIF;
E1_full = params.E1_full;
accel_factor = params.accel_factor;
seq_type = params.seq_type;
LUTCorrection = params.LUTCorrection;

AIF_s1 = perfusion_AIF_signal(squeeze(AIF(:,:,1,:)), mask(:,:,end), numPD, percentageAIF);
AIF_s2 = perfusion_AIF_signal(squeeze(AIF(:,:,2,:)), mask(:,:,end), numPD, percentageAIF);

AIF_s1_Full = perfusion_AIF_signal(squeeze(AIF_Full(:,:,1,:)), mask_Full(:,:,end), numPD, percentageAIF);
AIF_s2_Full = perfusion_AIF_signal(squeeze(AIF_Full(:,:,2,:)), mask_Full(:,:,end), numPD, percentageAIF);

figure; hold on; plot(AIF_s1); plot(AIF_s2, 'r'); plot(AIF_s1_Full, 'b-.'); plot(AIF_s2_Full, 'r-.'); hold off
legend('AIF mini, echo 0', 'AIF mini, echo 1', 'AIF full, echo 0', 'AIF full, echo 1');
title('Origianl AIF mini/full');

% t2* correction
[R, r, R2] = PerfusionT2StarCorrection_L1(AIF_s1(numPD+1:end), AIF_s2(numPD+1:end), TE0, TE1, 0.005);
R(1) = 1;
AIF_s1_R2Star = R .* AIF_s1(numPD+1:end);

[R_Full, r_Full, R2_Full] = PerfusionT2StarCorrection_L1(AIF_s1_Full(numPD+1:end), AIF_s2_Full(numPD+1:end), TE0, TE1);
R_Full(1) = 1;
AIF_s1_Full_R2Star = R_Full .* AIF_s1_Full(numPD+1:end);

figure; hold on; plot(AIF_s1(numPD+1:end)); plot(AIF_s1_R2Star); plot(AIF_s1_Full(numPD+1:end), 'b-.'); plot(AIF_s1_Full_R2Star, 'b-.'); hold off
legend('AIF mini, echo 0', 'AIF mini, echo 0, R2Star correction', 'AIF full, echo 0', 'AIF full, echo 0 , R2Star correction');
title('AIF mini/full with R2Star correction');

% t1 correction
if(LUTCorrection)
    Gd = [0:0.05:20]';
    [LUT, SR_Gd] = Matlab_gt_perfusion_bloch_simulation(single(Gd), FA_PD, FA_SR, TD, TR, E1_full, accel_factor, seq_type, seq_type, single(AIF_s1_R2Star(:)/AIF_s1(1)));
%     SR_Gd = interp1(LUT, Gd, AIF_s1_R2Star(:)/AIF_s1(1), 'linear');
    [LUT_Full, SR_Gd_Full] = Matlab_gt_perfusion_bloch_simulation(single(Gd), FA_PD, FA_SR, TD, TR, E1_full, accel_factor, seq_type, seq_type, single(AIF_s1_Full_R2Star(:)/AIF_s1(1)));
else
    ratio = sin(FA_PD*pi/180)/sin(FA_SR*pi/180);
    SR_Gd = AIF_s1_R2Star/AIF_s1(1) * ratio;
    SR_Gd = log(1-SR_Gd) * (-1/TS);
    
    SR_Gd_Full = AIF_s1_Full_R2Star/AIF_s1(1) * ratio;
    SR_Gd_Full = log(1-SR_Gd_Full) * (-1/TS);
    
end

figure; hold on; plot(SR_Gd); plot(SR_Gd_Full, 'b-.'); plot(SR_Gd*10, 'r'); hold off
legend('Gd mini', 'Gd full', '10*Gd mini');
title('Gd mini/full with R2Star and T1 correction');
