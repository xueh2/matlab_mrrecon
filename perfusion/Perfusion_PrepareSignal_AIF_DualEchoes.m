
function [cin, cin_noR2Star, cin_e1_noR2Star, cin_e0_R2Star, cin_e0, cin_e1, cin_PD] = Perfusion_PrepareSignal_AIF_DualEchoes(AIF_Echo0, AIF_Echo1, AIF_PD, mask, R2StarCorrection, TE0, TE1, FA_SR, FA_PD, TR, TS, TD, E1_full, accel_factor, seq_type, T1_0, T2_0, r1, r2, LUTCorrection)
% perform the aif cin function computation
% all times are in ms
% R2StarCorrection: whether to perform R2StarCorrection correction
% TE0, TE1: echo time for two echoes
% FA_SR, FA_PD : flip angle for AIF SR and PD images
% TR, TS, TD : aif readout TR, saturation time (from the end of SR pulse to the center kspace line) and TD (from end of SR pulse to the binning of first readout)
% E1_full, acce_factor
% seq_type : SSFP or Flash
% Gd_method : TwoPoint or LUT

%% compute mean PD

PD = squeeze(AIF_PD(:,:,:,:)); % first echo, first PD

RO = size(AIF_Echo0, 1);
E1 = size(AIF_Echo0, 2);

REP = size(AIF_Echo0, 3);

% get the raw cin
cin_e0 = zeros(REP, 1);
cin_e1 = zeros(REP, 1);

ind = find(mask(:)>0);

pd2d = PD(:,:,1);
vPD = pd2d(ind(:));

for r=1:REP
    im = AIF_Echo0(:,:,r);
    v = im(ind(:));
    v = v./vPD;
    cin_e0(r) = mean(v);
    
    im = AIF_Echo1(:,:,r);
    v = im(ind(:));
    v = v./vPD;
    cin_e1(r) = mean(v);
end

figure; hold on; plot(cin_e0); plot(cin_e1, 'r'); hold off

% get the PD
% for pp=1:size(PD, 3)
%     pd2d = PD(:,:,pp);
%     v = pd2d(ind(:));
%     cin_PD(pp) = mean(v(:));
% end
% cin_PD = max(cin_PD(:));

cin_PD = vPD;

% R2* correction
cin_e0_R2Star = cin_e0;
if(R2StarCorrection)
    [R, r, R2] = PerfusionT2StarCorrection_L1(cin_e0, cin_e1, TE0, TE1, 0.001);
    R(1) = 1;
    cin_e0_R2Star = R .* cin_e0;
end

figure; hold on; plot(cin_e0); plot(cin_e0_R2Star, 'b-.'); plot(cin_e1, 'r'); hold off
legend('echo 0', 'echo 0 with R2* correction', 'echo 1');

% T1 correction
if(LUTCorrection)
    Gd = [0:0.05:40]';
    % [LUT, SR_Gd] = Matlab_gt_perfusion_bloch_simulation(single(Gd), FA_PD, FA_SR, TD, TR, E1_full, accel_factor, seq_type, seq_type, T1_0, T2_0, single(cin_e0_R2Star(:)/cin_PD));
    [LUT, SR_Gd_R2Star] = Matlab_gt_perfusion_bloch_simulation(Gd, FA_PD, FA_SR, TD, TR, E1_full, accel_factor, seq_type, seq_type, T1_0, T2_0, r1, r2, single(cin_e0_R2Star(:)));
    
    [LUT, SR_Gd] = Matlab_gt_perfusion_bloch_simulation(Gd, FA_PD, FA_SR, TD, TR, E1_full, accel_factor, seq_type, seq_type, T1_0, T2_0, r1, r2, single(cin_e0(:)));
    
    [LUT, SR_Gd_e1] = Matlab_gt_perfusion_bloch_simulation(Gd, FA_PD, FA_SR, TD, TR, E1_full, accel_factor, seq_type, seq_type, T1_0, T2_0, r1, r2, single(cin_e1(:)));
else
    ratio = sin(FA_PD*pi/180)/sin(FA_SR*pi/180);
    SR_Gd = cin_e0_R2Star/cin_PD(1) * ratio;
    SR_Gd = log(1-SR_Gd) * (-1/TS);
      
end

cin = SR_Gd_R2Star;
cin_noR2Star = SR_Gd;
cin_e1_noR2Star = SR_Gd_e1;

