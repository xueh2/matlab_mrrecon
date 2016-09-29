function [LUT, Gd, SRcenter, PDcenter, SR, PD] = perfusion_lut_flash_pd_flash_sr(params, Gd);
% function [LUT, Gd] = perfusion_lut_flash_pd_flash_sr(params);
%
% calculation of look-up-table (LUT) for linearizing Magnetization (Mxy) vs [Gd]
%    PD flash
%    SR ssfp
% inputs:
%     GdDTPA relaxivities, r1, r2
%     initial T1, T2 without Gd (i.e., [Gd]=0)
%     sequence parameters (PD FA, SR FA, TR, TD, number of PE lines to center
%     of kspace)


% params.TR = 0.0027;
    % params.TR = 0.00112; % 1.12 ms (1-echo AIF protocol)
    % params.TR = 0.00212; % 2.12 ms (2 echo AIF protocol)
% params.TD = 0.005;
% params.PD_flip_angle = 8;
% params.SR_flip_angle = 50; 
% params.offresonance = 0;
% params.Npe_full = 111;
% params.PAT_accel_factor = 3;
% params.steadystate_prep = 'linear';
% params.N_runup = 5;
% params.rf_phase_spoiler_increment = 112;

% Gd = params.Gd; % Gd concentration (mM)

if(nargin < 2)
    Gd = [0:.1:20]; % Gd concentration (mM)
end

r1 = params.r1;
r2 = params.r2;
% T1_0 = 2.0; % s
% T2_0 = 0.250; % s
T1_0 = params.T1_0; % s
T2_0 = params.T2_0; % s

Npe_total = floor(params.Npe_full/params.PAT_accel_factor);
% Npe_center = round(Npe_total/2); % # PE lines to center of k-space
Npe_center = floor(Npe_total/2); % # PE lines to center of k-space
params.Npulses = Npe_center;


% PD signal (at [Gd]=0 without SR)
params.flip_angle = params.PD_flip_angle;
params.T1 = Gd2T1(0, T1_0, r1);
params.T2 = Gd2T1(0, T2_0, r2);
[PD] = flash_readout(params);
PDmag = squeeze(rss(PD(1:2,:,:),1)); % PE, GD
PDcenter = PDmag(end);

% SR signal vs [Gd]
params.flip_angle = params.SR_flip_angle;
M = zeros(4,params.Npulses,length(Gd));
for k = 1:length(Gd)
    params.T1 = Gd2T1(Gd(k), T1_0, r1);
    if(params.T1<0) params.T1 = abs(params.T1)+T1_0; end
    params.T2 = Gd2T1(Gd(k), T2_0, r2);
    if(params.T2<0) params.T2 = abs(params.T2)+T2_0; end
    M0 = 1 - exp(-params.TD/params.T1); % Saturation recovery initial magnetization at trigger delay, TD
    params.M0 = [0 0 M0 1]';
    [M(:,:,k)] = flash_readout(params);
end

% plot of Magn vs [Gd] for various TD at a fixed FA.
Mxy = squeeze(rss(M(1:2,:,:),1)); % PE, GD 
SRcenter = Mxy(end,:);
% figure; plot(Gd,SRcenter); hold on

LUT = SRcenter/PDcenter;

if nargout >= 5
    SR = SRcenter;
end
if nargout == 6
    PD = PDcenter;
end

% figure; plot(Gd, LUT)
% g = linspace(0,max(Gd),1000);
% [p,S] = polyfit(Gd(1:3),LUT(1:3),1);
% LUT_fit = p(1)*g + p(2);
% hold on;
% plot(g, LUT_fit,'m')
% simple_calculation = (sin(params.SR_flip_angle*pi/180)/sin(params.PD_flip_angle*pi/180))*(1 - exp(-(params.TD + Npe_center*params.TR)./Gd2T1(g, T1_0, r1)));
% plot(g,simple_calculation,':')
% grid on
% 

return







