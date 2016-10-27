
function [aif_TS, aif_TD, Gd_aif_slice, Gd_aif_slice_noR2Star, TS, TD, Gd_perf_slice] = PerformGadgetronRecon_GdPhantom_Test(h5Name, res_dir, debug_dir, r1, r2, T1_0_blood, T2_0_blood, T1_0_myo, T2_0_myo, startTube, Gd_tubes, perf_roi_file, Gd_aif_tubes, aif_roi_file, contrast, sliceprofile, sliceprofile_aif, B1)
% PerformGadgetronRecon_GdPhantom_Test(h5Name, res_dir, debug_dir, r1, r2, T1_0_blood, T2_0_blood, T1_0_myo, T2_0_myo, startTube, Gd_tubes, perf_roi_file, Gd_aif_tubes, aif_roi_file, contrast, sliceprofile, sliceprofile_aif, B1)

closeall
%% get perf parameters

aif_echo_time_0 = 0.65;
aif_echo_time_1 = 1.65;
aif_FA_PD = 8;
aif_FA_Perf = 8;
aif_TS = 21.17;
aif_TD = 4.21;
aif_TR = 2.12;
aif_E1_full = 35;
aif_accel_factor = 2;
aif_N_runup = 3;
aif_seq_type = 'Flash';
aif_Gd_method = 'LUT';
aif_PD_mean = 0;

FA_PD = 5;
FA_Perf = 49.9999;
TS = 65;
TD = 32.5;
TR = 2.5;
E1_full = 57;
accel_factor = 3;
N_runup = 3;
seq_type = 'SSFP';
seq_type_PD = 'Flash';
Gd_method = 'LUT';
% T1_0_blood = 1650;
% T2_0_blood = 220;
% T1_0_myo = 1175;
% T2_0_myo = 45;
post_T1_0_blood = 600;
post_T1_0_myo = 700;
check_contrast_status = 1;

dset = ismrmrd.Dataset(h5Name, 'dataset');
header = ismrmrd.xml.deserialize(dset.readxml());
dset.close();

FA_Perf = header.sequenceParameters.flipAngle_deg(1);
NumOfPD = header.userParameters.userParameterLong(2).value;

if(numel(header.userParameters.userParameterLong)==2)
    NumOfConcatenations = 1;
else
    NumOfConcatenations = header.userParameters.userParameterLong(3).value;
end

seq_type = header.sequenceParameters.sequence_type
header.acquisitionSystemInformation.systemFieldStrength_T

if(header.acquisitionSystemInformation.systemFieldStrength_T>2)
    if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
        FA_Perf = 40 * ((header.userParameters.userParameterDouble(6).value / header.userParameters.userParameterDouble(5).value) / 0.8436)
    else
        FA_Perf = 14 * ((header.userParameters.userParameterDouble(6).value / header.userParameters.userParameterDouble(5).value) / 0.44433)
    end
else
    if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
        FA_Perf = 50 * ((header.userParameters.userParameterDouble(6).value / header.userParameters.userParameterDouble(5).value) / (211.587 / 200))
    else
        FA_Perf = 14 * ((header.userParameters.userParameterDouble(6).value / header.userParameters.userParameterDouble(5).value) / (88.866 / 200))
    end
end

if(FA_Perf>header.sequenceParameters.flipAngle_deg(1))
    FA_Perf = header.sequenceParameters.flipAngle_deg(1);
end

aif_echo_time_0 = header.sequenceParameters.TE(2);
aif_echo_time_1 = header.sequenceParameters.TE(3);
aif_FA_PD = header.sequenceParameters.flipAngle_deg(2);
aif_FA_Perf = header.sequenceParameters.flipAngle_deg(2);
if(header.sequenceParameters.TI(2)==0)
    aif_TR = 2.12;
else
    aif_TR = header.sequenceParameters.echo_spacing(2);
end
% aif_E1_full = header.encoding(2).encodingLimits.kspace_encoding_step_1.maximum;
aif_E1_full = 2 * header.encoding(2).encodingLimits.kspace_encoding_step_1.center;
aif_accel_factor = header.encoding(2).parallelImaging.accelerationFactor.kspace_encoding_step_1;

while(mod(aif_E1_full, aif_accel_factor)~=0)
    aif_E1_full = aif_E1_full + 1;
end

if(header.sequenceParameters.TI(2)==0)
    aif_TD = 4.7;
    aif_TS = aif_TD + floor( floor(aif_E1_full / aif_accel_factor) / 2.0 ) * aif_TR;
else
    aif_TS = header.sequenceParameters.TI(2);
    aif_TD = aif_TS - floor( floor(aif_E1_full / aif_accel_factor) / 2.0 ) * aif_TR;
end

FA_PD = 5;
TR = header.sequenceParameters.echo_spacing(1);
% E1_full = header.encoding(1).encodedSpace.matrixSize.y;
% E1_full = 2*(header.encoding(1).encodingLimits.kspace_encoding_step_1.maximum - header.encoding(1).encodingLimits.kspace_encoding_step_1.center) + 1;
E1_full = 2*header.encoding(1).encodingLimits.kspace_encoding_step_1.center + 1;
accel_factor = header.encoding(1).parallelImaging.accelerationFactor.kspace_encoding_step_1;
while(mod(E1_full, accel_factor)~=0)
    E1_full = E1_full + 1;
end

if(header.sequenceParameters.TI(2)==0)
    if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
        TD = 17;
        TS = TD + TR * floor( floor(E1_full / accel_factor) / 2.0) - N_runup * TR;
    else
        TD = 8.5;
        TS = TD + TR * floor( floor(E1_full / accel_factor) / 2.0);
    end
else
    TS = header.sequenceParameters.TI(1);
    if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
        TD = TS - TR * floor( floor(E1_full / accel_factor) / 2.0) - N_runup * TR;
    else
        TD = TS - TR * floor( floor(E1_full / accel_factor) / 2.0);
    end
end

SLC = header.encoding(1).encodingLimits.slice.maximum+1;

disp(['--> Sequence parameters, AIF ... ']);
disp(['aif_echo_time_0 = ' num2str(aif_echo_time_0)]);
disp(['aif_echo_time_1 = ' num2str(aif_echo_time_1)]);
disp(['aif_FA_PD = ' num2str(aif_FA_PD)]);
disp(['aif_FA_Perf = ' num2str(aif_FA_Perf)]);
disp(['aif_TS = ' num2str(aif_TS)]);
disp(['aif_TD = ' num2str(aif_TD)]);
disp(['aif_TR = ' num2str(aif_TR)]);
disp(['aif_E1_full = ' num2str(aif_E1_full)]);
disp(['aif_accel_factor = ' num2str(aif_accel_factor)]);
disp(['aif_N_runup = ' num2str(aif_N_runup)]);
disp(['aif_seq_type = ' num2str(aif_seq_type)]);
disp(['aif_Gd_method = ' num2str(aif_Gd_method)]);
disp(['aif_PD_mean = ' num2str(aif_PD_mean)]);
disp(['%==================================================%']);
disp(['--> Sequence parameters, SR ... ']);
disp(['FA_PD = ' num2str(FA_PD)]);
disp(['FA_Perf = ' num2str(FA_Perf)]);
disp(['TS = ' num2str(TS)]);
disp(['TD = ' num2str(TD)]);
disp(['TR = ' num2str(TR)]);
disp(['E1_full = ' num2str(E1_full)]);
disp(['accel_factor = ' num2str(accel_factor)]);
disp(['N_runup = ' num2str(N_runup)]);
disp(['seq_type = ' num2str(seq_type)]);
disp(['seq_type_PD = ' num2str(seq_type_PD)]);
disp(['Gd_method = ' num2str(Gd_method)]);
disp(['T1_0_blood = ' num2str(T1_0_blood)]);
disp(['T2_0_blood = ' num2str(T2_0_blood)]);
disp(['T1_0_myo = ' num2str(T1_0_myo)]);
disp(['T2_0_myo = ' num2str(T2_0_myo)]);
disp(['post_T1_0_blood = ' num2str(post_T1_0_blood)]);
disp(['post_T1_0_myo = ' num2str(post_T1_0_myo)]);
disp(['check_contrast_status = ' num2str(check_contrast_status)]);
disp(['%==================================================%']);

%% load data
Im_TI = readGTPlusExportImageSeries_Squeeze(res_dir, 103);
cd(debug_dir)

SR_Norm_TI = analyze75read('SRNorm_0.hdr');
Gd_TI = analyze75read('Input_perf_computeFlowMap_0.hdr');

if(strcmp(lower(seq_type), 'ssfp'))
    LUT_TI = analyze75read('Perf_T1_Correction_LUT_perf_ssfp_PD_flash.hdr');
else
    LUT_TI = analyze75read('Perf_T1_Correction_LUT_perf_flash_PD_flash.hdr');
end

LUT_aif_TI = analyze75read('aif_cin_LUT_flash_pd_flash_sr.hdr');
aif_0 = analyze75read('input_aif_0.hdr'); size(aif_0)
aif_1 = analyze75read('input_aif_1.hdr');
aif_pd_0 = analyze75read('input_aif_PD_0.hdr'); size(aif_pd_0)
aif_pd_1 = analyze75read('input_aif_PD_1.hdr');

save AIF aif_0 aif_1 aif_pd_0 aif_pd_1 Im_TI

%% convert to Gd

figure; imagescn(Im_TI, [], [], [], 3);

% m = zeros(18, 24);

mask_file = perf_roi_file;
% for ii=1:24  
%     v = roi_timeseries(Im_TI, mask_file, ii, 1);
%     p = v.m;
%     m(:, ii) = p(:);      
% end

SR_PD = zeros(3, 24);
Im = Im_TI(:,:,8:end);
Im = mean(Im, 3);
srpd = Im ./ Im_TI(:,:,1);
for ii=1:24  
    v = roi_timeseries(srpd, mask_file, ii, 1);
    p = v.m;
    SR_PD(1, ii) = p(:);      
end

% plot LUT

Gd_LUT = [0:0.01:5];
[LUT_TI_m, signal_SR, signal_PD] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT', FA_PD, FA_Perf, TD, TR, E1_full, accel_factor, seq_type, seq_type_PD, T1_0_myo, T2_0_myo, r1, r2, []);

params.r1 = r1;
params.r2 = r2;
params.T1_0 = T1_0_myo/1e3;
params.T2_0 = T2_0_myo/1e3;
params.TR = TR/1e3;
params.TD = TD/1e3;
params.PD_flip_angle = FA_PD;
params.SR_flip_angle = FA_Perf;
params.offresonance = 0;
params.Npe_full = E1_full;
params.PAT_accel_factor =  accel_factor;
params.steadystate_prep = 'linear';
params.N_runup = 3;
params.rf_phase_spoiler_increment = 112;

% compute Gd for every tube, taking into account of B1 map and slice profile

perf_Gd = zeros(1, 24);
perf_Gd_slice = zeros(1, 24);

if(isempty(B1))
    
    LUTslice = [];
    if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
        [LUTslice, Gd] = perfusion_lut_flash_pd_ssfp_sr(params, Gd_LUT);
        LUT_TI_mm = LUTslice;
    else   
        % [LUT_TI_mm, Gd] = perfusion_lut_flash_pd_flash_sr(params, Gd_LUT);

        [LUT_TI_mm, SR, PD] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT', FA_PD, FA_Perf, TD, TR, E1_full, accel_factor, seq_type, seq_type_PD, T1_0_myo, T2_0_myo, r1, r2, []);

        if(~isempty(sliceprofile))
            clear LUT SRc PDc SR PD
            for i = 1:length(sliceprofile)
                disp(['Slice profile : ' num2str(i)])
                params.PD_flip_angle = FA_PD * sliceprofile(i);
                params.SR_flip_angle = FA_Perf * sliceprofile(i);
                [LUT(i,:), SR(i,:), PD(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT', params.PD_flip_angle, params.SR_flip_angle, TD, TR, E1_full, accel_factor, seq_type, seq_type_PD, T1_0_myo, T2_0_myo, r1, r2, []);
            end

            LUTslice = sum(cat(1,SR,SR(2:end,:)),1)/sum(cat(1,PD,PD(2:end,:)),1);
        end
    end
    
    for tt=startTube:startTube+numel(Gd_tubes)-1
        sr_over_pd = SR_PD(1, tt);
        perf_Gd_slice(1, tt) = interp1(LUTslice, Gd_LUT, sr_over_pd);  
        perf_Gd(1, tt) = interp1(LUT_TI_m, Gd_LUT, sr_over_pd);
    end 
    
    figure;
    hold on

    plot(Gd_LUT, LUT_TI_m, 'b--');

    if(~isempty(LUTslice))
        plot(Gd_LUT, LUTslice, 'k');
    end

    plot(Gd_LUT, LUT_TI_mm, 'g-.');

    hold off
    if(~isempty(LUTslice))
        legend(['TD=' num2str(TD)], ['with slice profile']);
    else
        legend(['TD=' num2str(TD)]);
    end
    title('Perf')
    xlabel('Gd')
    ylabel('SR/PD')
else
    for tt=startTube:startTube+numel(Gd_tubes)-1

        title_str = ['Perf - processing tube : ' num2str(tt) ' - B1 ' num2str(B1(tt))];
        disp(title_str);

        [LUT_TI_mm, SR, PD] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT', FA_PD, FA_Perf, TD, TR, E1_full, accel_factor, seq_type, seq_type_PD, T1_0_myo, T2_0_myo, r1, r2, []);
        
        LUTslice = [];
        if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
            params.PD_flip_angle = FA_PD * B1(tt);
            params.SR_flip_angle = FA_Perf * B1(tt);           
            [LUTslice, Gd] = perfusion_lut_flash_pd_ssfp_sr(params, Gd_LUT);
        else   
            % [LUT_TI_mm, Gd] = perfusion_lut_flash_pd_flash_sr(params, Gd_LUT);

            % [LUT_TI_mm, SR, PD] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT', FA_PD, FA_Perf, TD, TR, E1_full, accel_factor, seq_type, seq_type_PD, T1_0_myo, T2_0_myo, r1, r2, []);

            if(~isempty(sliceprofile))
                clear LUT SRc PDc SR PD
                for i = 1:length(sliceprofile)
                    disp(['Slice profile : ' num2str(i)])
                    params.PD_flip_angle = FA_PD * sliceprofile(i) * B1(tt);
                    params.SR_flip_angle = FA_Perf * sliceprofile(i) * B1(tt);
                    % [LUT(i,:), Gd, SRc, PDc, SR(i,:), PD(i,:)] = perfusion_lut_flash_pd_flash_sr(params, Gd_LUT);

                    [LUT(i,:), SR(i,:), PD(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT', params.PD_flip_angle, params.SR_flip_angle, TD, TR, E1_full, accel_factor, seq_type, seq_type_PD, T1_0_myo, T2_0_myo, r1, r2, []);
                end

                LUTslice = sum(cat(1,SR,SR(2:end,:)),1)/sum(cat(1,PD,PD(2:end,:)),1);
            end
        end

        sr_over_pd = SR_PD(1, tt);
        perf_Gd_slice(1, tt) = interp1(LUTslice, Gd_LUT, sr_over_pd);  
        perf_Gd(1, tt) = interp1(LUT_TI_m, Gd_LUT, sr_over_pd);   

        % plot
        figure;
        hold on

        plot(Gd_LUT, LUT_TI_m, 'b--');

        if(~isempty(LUTslice))
            plot(Gd_LUT, LUTslice, 'k');
        end

        plot(Gd_LUT, LUT_TI_mm, 'g-.');

        hold off
        if(~isempty(LUTslice))
            legend(['TD=' num2str(TD)], ['with slice profile']);
        else
            legend(['TD=' num2str(TD)]);
        end
        title(title_str)
        xlabel('Gd')
        ylabel('SR/PD')
    end
end

% % compute Gd
% perf_Gd = zeros(1, 24);
% for ii=1:size(perf_Gd, 2)       
%     sr_over_pd = SR_PD(1, ii);
%     perf_Gd(1, ii) = interp1(LUT_TI_m, Gd_LUT, sr_over_pd);   
% end
% 
% if(~isempty(LUTslice))
%     perf_Gd_slice = zeros(1, 24);
%     for ii=1:size(perf_Gd_slice, 2)       
%         sr_over_pd = SR_PD(1, ii);
%         perf_Gd_slice(1, ii) = interp1(LUTslice, Gd_LUT, sr_over_pd);   
%     end
% end

Gd = Gd_tubes;

Gd_perf = perf_Gd(:, startTube:startTube+numel(Gd)-1);
P = polyfit(Gd(:), Gd_perf(:),1);
Gd_fit = P(1)*Gd+P(2);

xpos = Gd(2)
ypos = 0.85 * max(Gd_fit)
theString = sprintf('y = %.5f x + %.5f', P(1), P(2));

if(~isempty(LUTslice))
    Gd_perf_slice = perf_Gd_slice(:, startTube:startTube+numel(Gd)-1);
    P = polyfit(Gd(:), Gd_perf_slice(:),1);
    Gd_fit_slice = P(1)*Gd+P(2);

    xpos_slice = Gd(2)
    ypos_slice = 0.65 * max(Gd_fit_slice)
    theString_slice = sprintf('y = %.5f x + %.5f', P(1), P(2));
end

figure; 
hold on
plot(Gd, Gd_perf(1,:), 'bx'); 

if(~isempty(LUTslice))
    plot(Gd, Gd_perf_slice(1,:), 'k+');
    plot(Gd, Gd_fit_slice(:), 'k--'); 
    text(xpos_slice, ypos_slice, theString_slice, 'FontSize', 12);
end

plot(Gd, Gd_fit(:), 'r-');
text(xpos, ypos, theString, 'FontSize', 12);

hold off
if(~isempty(LUTslice))
    legend(['TD=' num2str(TD)], 'with slice profile');
else
    legend(['TD=' num2str(TD)]);
end

title([ contrast ', perfusion, Gd conversion, r1=' num2str(r1) ' - r2=' num2str(r2)])
xlabel('Gd, truth'); 
ylabel('Gd, perf');
box on
grid on

% plot 2
figure; 
hold on
plot(Gd, Gd_perf_slice(1,:), 'k+');
plot(Gd, Gd_fit_slice(:), 'k--'); 
text(xpos_slice, ypos_slice, theString_slice, 'FontSize', 12);

hold off
legend(['TD=' num2str(TD)]);

title([ contrast ', perfusion, Gd conversion, r1=' num2str(r1) ' - r2=' num2str(r2)])
xlabel('Gd, truth'); 
ylabel('Gd, perf');
box on
grid on

%% aif

figure; imagescn(aif_pd_0(:,:,1));

figure; imagescn(cat(4, aif_pd_0, aif_pd_1), [], [], [], 3);

figure; imagescn(cat(4, aif_0, aif_1), [], [], [], 3);

N = size(aif_0, 3);
Npd = size(aif_pd_0, 3);

mask_file = aif_roi_file;

aif_echo0 = aif_0;
aif_echo1 = aif_1;

aif_echo0_PD = aif_pd_0;
aif_echo0_PD = aif_echo0_PD(:,:,1);

aif_echo1_PD = aif_pd_1;
aif_echo1_PD = aif_echo1_PD(:,:,1);

aif_echo0 = mean(aif_echo0(:,:,7:end), 3);
aif_echo1 = mean(aif_echo1(:,:,7:end), 3);

I = cat(3, aif_echo0, aif_echo1);

figure; imagescn(I, [], [], [], 3); 

params.r1 = r1;
params.r2 = r2;
params.T1_0 = T1_0_blood/1e3;
params.T2_0 = T2_0_blood/1e3;
params.TR = aif_TR/1e3;
params.TD = aif_TD/1e3;
params.PD_flip_angle = aif_FA_PD;
params.SR_flip_angle = aif_FA_Perf;
params.offresonance = 0;
params.Npe_full = aif_E1_full;
params.PAT_accel_factor =  aif_accel_factor;
params.steadystate_prep = 'linear';
params.N_runup = 3;
params.rf_phase_spoiler_increment = 112;
    
Gd_LUT = [0:0.01:14]';
[LUT, signal_SR, signal_PD] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, aif_FA_PD, aif_FA_Perf, aif_TD, aif_TR, aif_E1_full, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);
% [LUT, Gd, SRcenter, PDcenter, SR, PD] = perfusion_lut_flash_pd_flash_sr(params, Gd_LUT);

aif_Gd = zeros(7, 24);
R2Star = zeros(1, 24);
cin_e0 = zeros(1, 24);
cin_e1 = zeros(1, 24);
PDv = zeros(1, 24);
cin_e0_R2Star = zeros(1, 24);  

if(isempty(B1))
    params.r1 = r1;
    params.r2 = r2;
    params.T1_0 = T1_0_blood/1e3;
    params.T2_0 = T2_0_blood/1e3;
    params.TR = aif_TR/1e3;
    params.TD = aif_TD/1e3;
    params.PD_flip_angle = aif_FA_PD;
    params.SR_flip_angle = aif_FA_Perf;
    params.offresonance = 0;
    params.Npe_full = aif_E1_full;
    params.PAT_accel_factor =  aif_accel_factor;
    params.steadystate_prep = 'linear';
    params.N_runup = 3;
    params.rf_phase_spoiler_increment = 112;

    clear LUTOne SRc PDc SR PD

    SR = zeros(length(sliceprofile_aif), numel(Gd_LUT));
    SR1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
    SR2 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
    SR3 = zeros(length(sliceprofile_aif), numel(Gd_LUT));

    PD = zeros(length(sliceprofile_aif), numel(Gd_LUT));
    PD1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
    PD2 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
    PD3 = zeros(length(sliceprofile_aif), numel(Gd_LUT));

    % without B1
    SR_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
    SR1_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
    SR2_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));

    PD_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
    PD1_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
    PD2_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));

    for i = 1:length(sliceprofile_aif)            
        params.PD_flip_angle = aif_FA_PD * sliceprofile_aif(i);
        params.SR_flip_angle = aif_FA_Perf * sliceprofile_aif(i);

        disp(['AIF slice profile : ' num2str(i), ' - FA_PD ' num2str(params.PD_flip_angle) ' - FA_SR ' num2str(params.SR_flip_angle)])

        % [LUTOne(i,:), Gd, SRc, PDc, SR(i,:), PD(i,:)] = perfusion_lut_flash_pd_flash_sr(params, Gd_LUT);
        % [LUTOne(i,:), SR(i,:), PD(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);

        [LUTOne(i,:), SR(i,:), PD(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);
        [LUTOne(i,:), SR1(i,:), PD1(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full+2*aif_accel_factor, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);
        [LUTOne(i,:), SR2(i,:), PD2(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full-2*aif_accel_factor, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);
        [LUTOne(i,:), SR3(i,:), PD3(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full+4*aif_accel_factor, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);
    end

    % with B1, slice profile, averaging
%         SR_all = SR + SR1 + SR2;
%         PD_all = PD + PD1 + PD2;       

    SR_all = SR + SR1;
    PD_all = PD + PD1;
    LUT_slice_averaging = sum(cat(1,SR_all,SR_all(2:end,:)),1)./sum(cat(1,PD_all,PD_all(2:end,:)),1);

    % with slice profile only
    SR_all = SR;
    PD_all = PD;
    LUT_slice = sum(cat(1,SR_all,SR_all(2:end,:)),1)./sum(cat(1,PD_all,PD_all(2:end,:)),1);

    % with averaging only
    SR_all = SR(1,:) + SR1(1,:);
    PD_all = PD(1,:) + PD1(1,:);
    LUT_averaging = SR_all./PD_all;
    
    % plot
    figure;
    hold on
    plot(Gd_LUT, LUT, 'b--');

    plot(Gd_LUT, LUT_slice, 'k--');
    plot(Gd_LUT, LUT_averaging, 'm-.');
    plot(Gd_LUT, LUT_slice_averaging, 'Color', [0.5 0.5 0.5]);

    hold off
    legend('no slice, no B1, no averaging', 'only slice profile', 'only averaging', 'Slice and averaging');
    title('AIF')
    xlabel('Gd')
    ylabel('SR/PD')
    box on
    grid on
        
    for tt=startTube:startTube+numel(Gd_aif_tubes)-1 
        aif_Gd(1, tt) = tt;
        p1 = roi_timeseries(I, mask_file, tt, 1);
        
        pd1 = roi_timeseries(aif_echo0_PD(:,:,1), mask_file, tt, 1);
        pd2 = roi_timeseries(aif_echo1_PD(:,:,1), mask_file, tt, 1);

        cin_e0(tt) = p1.m(1);
        cin_e1(tt) = p1.m(2);        

        R2Star(tt) = log(pd1.m(1)/pd2.m(1)) / (aif_echo_time_1-aif_echo_time_0);
        if(R2Star(tt)<=0)
            R2Star(tt)=0;
        end
        R = exp(aif_echo_time_0*R2Star(tt));
        cin_e0_R2Star(tt) = R .* cin_e0(tt);

        PDv(tt) = R * pd1.m(1);
        
        sr_over_pd = cin_e0_R2Star(tt)/PDv(tt);
        aif_Gd(2, tt) = interp1(LUT, Gd_LUT, sr_over_pd); % no B1, no slice, no averaging, with T2*

        aif_Gd(4, tt) = interp1(LUT_averaging, Gd_LUT, sr_over_pd); % with B1 and slice and averaging, with T2*
        aif_Gd(6, tt) = interp1(LUT_slice, Gd_LUT, sr_over_pd); % with slice profile only

        sr_over_pd = cin_e0(tt)/PDv(tt);
        aif_Gd(3, tt) = interp1(LUT, Gd_LUT, sr_over_pd); % no B1, no slice, no averaging, no T2*
        aif_Gd(5, tt) = interp1(LUT_averaging, Gd_LUT, sr_over_pd); % with B1 and slice and averaging, no T2*
    end
else
    for tt=startTube:startTube+numel(Gd_aif_tubes)-1 

        title_str = ['AIF - processing tube : ' num2str(tt) ' - B1 ' num2str(B1(tt)) ' - TD ' num2str(aif_TD)];
        disp(title_str)

        if(~isempty(sliceprofile_aif))

            params.r1 = r1;
            params.r2 = r2;
            params.T1_0 = T1_0_blood/1e3;
            params.T2_0 = T2_0_blood/1e3;
            params.TR = aif_TR/1e3;
            params.TD = aif_TD/1e3;
            params.PD_flip_angle = aif_FA_PD;
            params.SR_flip_angle = aif_FA_Perf;
            params.offresonance = 0;
            params.Npe_full = aif_E1_full;
            params.PAT_accel_factor =  aif_accel_factor;
            params.steadystate_prep = 'linear';
            params.N_runup = 3;
            params.rf_phase_spoiler_increment = 112;

            clear LUTOne SRc PDc SR PD

            SR = zeros(length(sliceprofile_aif), numel(Gd_LUT));
            SR1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
            SR2 = zeros(length(sliceprofile_aif), numel(Gd_LUT));

            PD = zeros(length(sliceprofile_aif), numel(Gd_LUT));
            PD1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
            PD2 = zeros(length(sliceprofile_aif), numel(Gd_LUT));

            % without B1
            SR_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
            SR1_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
            SR2_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));

            PD_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
            PD1_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));
            PD2_noB1 = zeros(length(sliceprofile_aif), numel(Gd_LUT));

            for i = 1:length(sliceprofile_aif)            
                params.PD_flip_angle = aif_FA_PD * sliceprofile_aif(i) * B1(tt);
                params.SR_flip_angle = aif_FA_Perf * sliceprofile_aif(i) * B1(tt);

                disp(['AIF slice profile : ' num2str(i), ' - FA_PD ' num2str(params.PD_flip_angle) ' - FA_SR ' num2str(params.SR_flip_angle)])

                % [LUTOne(i,:), Gd, SRc, PDc, SR(i,:), PD(i,:)] = perfusion_lut_flash_pd_flash_sr(params, Gd_LUT);
                % [LUTOne(i,:), SR(i,:), PD(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);

                [LUTOne(i,:), SR(i,:), PD(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);
                [LUTOne(i,:), SR1(i,:), PD1(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full+2*aif_accel_factor, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);
                [LUTOne(i,:), SR2(i,:), PD2(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full-2*aif_accel_factor, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);

            end

            for i = 1:length(sliceprofile_aif)            
                params.PD_flip_angle = aif_FA_PD * sliceprofile_aif(i);
                params.SR_flip_angle = aif_FA_Perf * sliceprofile_aif(i);

                disp(['AIF slice profile : ' num2str(i), ' - FA_PD ' num2str(params.PD_flip_angle) ' - FA_SR ' num2str(params.SR_flip_angle)])

                [LUTOne(i,:), SR_noB1(i,:), PD_noB1(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);
                [LUTOne(i,:), SR1_noB1(i,:), PD1_noB1(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full+2*aif_accel_factor, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);
                [LUTOne(i,:), SR2_noB1(i,:), PD2_noB1(i,:)] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT, params.PD_flip_angle, params.SR_flip_angle, aif_TD, aif_TR, aif_E1_full-2*aif_accel_factor, aif_accel_factor, aif_seq_type, aif_seq_type, T1_0_blood, T2_0_blood, r1, r2, []);

            end

            % with B1, slice profile, averaging
    %         SR_all = SR + SR1 + SR2;
    %         PD_all = PD + PD1 + PD2;       

            SR_all = SR + SR1;
            PD_all = PD + PD1;       
            LUT_slice_B1_averaging = sum(cat(1,SR_all,SR_all(2:end,:)),1)./sum(cat(1,PD_all,PD_all(2:end,:)),1);

            % with B1, averaging
            SR_all = SR(1,:) + SR1(1,:);
            PD_all = PD(1,:) + PD1(1,:);
            LUT_B1_averaging = SR_all./PD_all;

            % with slice and averaging
            SR_all = SR_noB1 + SR1_noB1;
            PD_all = PD_noB1 + PD1_noB1;
            LUT_slice_averaging = sum(cat(1,SR_all,SR_all(2:end,:)),1)./sum(cat(1,PD_all,PD_all(2:end,:)),1);

            % with B1 only
            SR_all = SR(1,:);
            PD_all = PD(1,:);
            LUT_B1 = SR_all./PD_all;

            % with slice profile only
            SR_all = SR_noB1(1,:);
            PD_all = PD_noB1(1,:);
            LUT_slice = SR_all./PD_all;

            % with averaging only
            SR_all = SR_noB1(1,:) + SR1_noB1(1,:);
            PD_all = PD_noB1(1,:) + PD1_noB1(1,:);
            LUT_averaging = SR_all./PD_all;
        end

        % plot
        figure;
        hold on
        plot(Gd_LUT, LUT, 'b--');

        plot(Gd_LUT, LUT_slice_B1_averaging, 'r--');
        plot(Gd_LUT, LUT_slice, 'k--');
        plot(Gd_LUT, LUT_B1, 'g--');
        plot(Gd_LUT, LUT_averaging, 'm-');
        plot(Gd_LUT, LUT_B1_averaging, 'c-.');
        plot(Gd_LUT, LUT_slice_averaging, 'Color', [0.5 0.5 0.5]);

        hold off
        legend(['no slice, no B1, no averaging'], ['slice/B1/averaging'], 'only slice profile', 'only B1', 'only averaging', 'B1 and averaging', 'Slice and averaging');
        title(title_str)
        xlabel('Gd')
        ylabel('SR/PD')
        box on
        grid on

        % R2star

        aif_Gd(1, tt) = tt;
        p1 = roi_timeseries(I, mask_file, tt, 1);
        PD = roi_timeseries(aif_echo0_PD(:,:,1), mask_file, tt, 1);

        cin_e0(tt) = p1.m(1);
        cin_e1(tt) = p1.m(2);
        PDv(tt) = PD.m(1);

        R2Star(tt) = log(cin_e0(tt)/cin_e1(tt)) / (aif_echo_time_1-aif_echo_time_0);
        if(R2Star(tt)<=0)
            R2Star(tt)=0;
        end
        R = exp(aif_echo_time_0*R2Star(tt));
        cin_e0_R2Star(tt) = R .* cin_e0(tt);

        sr_over_pd = cin_e0_R2Star(tt)/PDv(tt);
        aif_Gd(2, tt) = interp1(LUT, Gd_LUT, sr_over_pd); % no B1, no slice, no averaging, with T2*

        aif_Gd(4, tt) = interp1(LUT_slice_B1_averaging, Gd_LUT, sr_over_pd); % with B1 and slice and averaging, with T2*
        aif_Gd(6, tt) = interp1(LUT_slice, Gd_LUT, sr_over_pd); % with slice profile only
        aif_Gd(7, tt) = interp1(LUT_B1, Gd_LUT, sr_over_pd); % with B1 only
        aif_Gd(8, tt) = interp1(LUT_averaging, Gd_LUT, sr_over_pd); % with averaging only

        sr_over_pd = cin_e0(tt)/PDv(tt);
        aif_Gd(3, tt) = interp1(LUT, Gd_LUT, sr_over_pd); % no B1, no slice, no averaging, no T2*
        aif_Gd(5, tt) = interp1(LUT_slice_B1_averaging, Gd_LUT, sr_over_pd); % with B1 and slice and averaging, no T2*
    end
end

Gd = Gd_aif_tubes;
Gd_aif = aif_Gd(2, startTube:startTube+numel(Gd)-1);
Gd_echo0_aif = aif_Gd(3, startTube:startTube+numel(Gd)-1);

if(~isempty(sliceprofile_aif))
    Gd_aif_slice = aif_Gd(4, startTube:startTube+numel(Gd)-1);
    
    P = polyfit(Gd(:), Gd_aif_slice(:),1);
    Gd_aif_slice_fit = P(1)*Gd+P(2);
    xpos_slice = Gd(2);
    ypos_slice = 0.45 * max(Gd_aif_slice_fit);
    theString_slice = sprintf('with T2* correction/B1Map/slice profile/averaging,\n y = %.5f x + %.5f', P(1), P(2));
    
    Gd_aif_slice_noR2Star = aif_Gd(5, startTube:startTube+numel(Gd)-1);
    P = polyfit(Gd(:), Gd_aif_slice_noR2Star(:),1);
    Gd_aif_slice_noR2Star_fit = P(1)*Gd+P(2);
    xpos_slice_noR2Star = Gd(2);
    ypos_slice_noR2Star = 0.45 * max(Gd_aif_slice_noR2Star_fit);
    theString_slice_noR2Star = sprintf('with B1Map/slice profile/averaging, no T2* correction,\n y = %.5f x + %.5f', P(1), P(2));
end

P = polyfit(Gd(:), Gd_aif(:),1);
Gd_fit = P(1)*Gd+P(2);
xpos = Gd(2);
ypos = 0.85 * max(Gd_fit);
theString = sprintf('with T2* correction, without B1Map/slice profile/averaging,\n y = %.5f x + %.5f', P(1), P(2));

P = polyfit(Gd(:), Gd_echo0_aif(:),1);
Gd_echo0_fit = P(1)*Gd+P(2);
xpos2 = Gd(2);
ypos2 = 0.65 * max(Gd_fit);
theString2 = sprintf('without T2* correction, without B1Map/slice profile/averaging,\n y = %.5f x + %.5f', P(1), P(2));

figure;
hold on
plot(Gd, Gd_aif, 'x'); 
plot(Gd, Gd_echo0_aif, 'r+'); 

if(~isempty(sliceprofile_aif))
    plot(Gd, Gd_aif_slice, 'o'); 
    plot(Gd, Gd_aif_slice_fit, 'k-');
    text(xpos_slice, ypos_slice, theString_slice, 'FontSize', 12);
    
    plot(Gd, Gd_aif_slice_noR2Star, 'o'); 
    plot(Gd, Gd_aif_slice_noR2Star_fit, 'g-');
    text(xpos_slice, ypos_slice_noR2Star, theString_slice_noR2Star, 'FontSize', 12);
end

plot(Gd, Gd_fit, 'b-'); 
plot(Gd, Gd_echo0_fit, 'r-.'); 
text(xpos, ypos, theString, 'FontSize', 12);
text(xpos2, ypos2, theString2, 'FontSize', 12);

hold off
if(~isempty(sliceprofile_aif))
    legend('With T2* correction, without B1Map/slice profile/averaging', 'Without T2* correction, without B1Map/slice profile/averaging', 'With T2* correction, with B1Map/slice profile/averaging', 'Without T2* correction, with B1Map/slice profile/averaging');
else   
    legend('With T2* correction', 'Without T2* correction');
end
title([contrast ', aif, Gd conversion, r1=' num2str(r1) ' - r2=' num2str(r2)])
xlabel('Gd, truth'); 
ylabel('Gd, aif');
box on
grid on

% fig 2
figure;
hold on

plot(Gd, Gd_aif_slice, 'o'); 
plot(Gd, Gd_aif_slice_fit, 'k-');
text(xpos_slice, ypos_slice, theString_slice, 'FontSize', 12);

plot(Gd, Gd_aif_slice_noR2Star, 'o'); 
plot(Gd, Gd_aif_slice_noR2Star_fit, 'g-');
text(xpos_slice, ypos_slice_noR2Star, theString_slice_noR2Star, 'FontSize', 12);
    
hold off
legend('With T2* correction', 'Without T2* correction');
title([contrast ', aif, Gd conversion, r1=' num2str(r1) ' - r2=' num2str(r2)])
xlabel('Gd, truth'); 
ylabel('Gd, aif');
box on
grid on


