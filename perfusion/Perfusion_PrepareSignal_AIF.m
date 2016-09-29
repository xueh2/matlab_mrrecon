
function [cin, cin_e0_R2Star, cin_PD] = Perfusion_PrepareSignal_AIF(AIF_Echo0, AIF_PD, mask, FA_SR, FA_PD, TR, TS, TD, E1_full, accel_factor, seq_type, LUTCorrection)
% perform the aif cin function computation
% all times are in ms
% R2StarCorrection: whether to perform R2StarCorrection correction
% FA_SR, FA_PD : flip angle for AIF SR and PD images
% TR, TS, TD : aif readout TR, saturation time (from the end of SR pulse to the center kspace line) and TD (from end of SR pulse to the binning of first readout)
% E1_full, acce_factor
% seq_type : SSFP or Flash
% Gd_method : TwoPoint or LUT

%% compute mean PD

PD = AIF_PD(:,:,1); % first echo, first PD

RO = size(AIF_Echo0, 1);
E1 = size(AIF_Echo0, 2);

REP = size(AIF_Echo0, 3);

% get the raw cin
cin_e0 = zeros(REP, 1);
cin_e1 = zeros(REP, 1);

ind = find(mask(:)>0);

for r=1:REP
    im = AIF_Echo0(:,:,r);
    v = im(ind(:));
    v = sort(v);
    cin_e0(r) = mean(v( ceil(0.75*numel(v)):end));   
end

figure; hold on; plot(cin_e0); hold off

% get the PD
v = PD(ind(:));
cin_PD = mean(v(:));

% R2* correction
cin_e0_R2Star = cin_e0;

figure; hold on; plot(cin_e0); plot(cin_e0_R2Star, 'b-.'); plot(cin_e1, 'r'); hold off
legend('echo 0', 'echo 0 with R2* correction', 'echo 1');

% T1 correction
if(LUTCorrection)
    Gd = [0:0.05:20]';
    [LUT, SR_Gd] = Matlab_gt_perfusion_bloch_simulation(single(Gd), FA_PD, FA_SR, TD, TR, E1_full, accel_factor, seq_type, seq_type, single(cin_e0_R2Star(:)/cin_PD));
else
    ratio = sin(FA_PD*pi/180)/sin(FA_SR*pi/180);
    SR_Gd = cin_e0_R2Star/cin_PD(1) * ratio;
    SR_Gd = log(1-SR_Gd) * (-1/TS);
      
end

cin = SR_Gd;
