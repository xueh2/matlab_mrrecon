
function [TS, TD, Gd_perf_slice] = PerformGadgetronRecon_T1Phantom_For_Gd_Test(h5Name, res_dir, debug_dir, r1, r2, T1_0_myo, T2_0_myo, startTube, Gd_tubes, perf_roi_file, sliceprofile)
% [TS, TD, Gd_perf_slice] = PerformGadgetronRecon_T1Phantom_For_Gd_Test(h5Name, res_dir, debug_dir, r1, r2, T1_0_myo, T2_0_myo, startTube, Gd_tubes, perf_roi_file, sliceprofile)

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
disp(['T1_0_myo = ' num2str(T1_0_myo(1))]);
disp(['T2_0_myo = ' num2str(T2_0_myo(1))]);
disp(['check_contrast_status = ' num2str(check_contrast_status)]);
disp(['%==================================================%']);

%% load data
Im_TI = readGTPlusExportImageSeries_Squeeze(res_dir, 103);
cd(debug_dir)

SR_Norm_TI = analyze75read('SRNorm_0.hdr');
Gd_TI = analyze75read('Input_perf_computeFlowMap_0.hdr');

if(strcmp(lower(seq_type), 'ssfp'))
    LUT_TI = analyze75read('Perf_T1_Correction_LUT_Valid.hdr');
else
    LUT_TI = analyze75read('Perf_T1_Correction_LUT_Valid.hdr');
end

LUT_aif_TI = analyze75read('aif_cin_LUT_Valid.hdr');
aif_0 = analyze75read('input_aif_0.hdr'); size(aif_0)
aif_1 = analyze75read('input_aif_1.hdr');
aif_pd_0 = analyze75read('input_aif_PD_0.hdr'); size(aif_pd_0)
aif_pd_1 = analyze75read('input_aif_PD_1.hdr');

save AIF aif_0 aif_1 aif_pd_0 aif_pd_1 Im_TI

%% convert to Gd

SLC = 1;
REP = size(Im_TI,3);
if(numel(size(Im_TI))==4)
    SLC = size(Im_TI, 3);
    REP = size(Im_TI, 4);
end

RO = size(Im_TI,1);
E1 = size(Im_TI,2);

Im_TI_all = reshape(Im_TI, [RO E1 SLC REP]);
Gd_perf_slice_all = [];

for slc=1:SLC

    Im_TI = squeeze(Im_TI_all(:,:,slc,:));
    
    figure; imagescn(Im_TI, [], [], [], 3);

    % m = zeros(18, 24);

    mask_file = perf_roi_file;
    % for ii=1:24  
    %     v = roi_timeseries(Im_TI, mask_file, ii, 1);
    %     p = v.m;
    %     m(:, ii) = p(:);      
    % end

    SR_PD = zeros(1, 2);
    Im = Im_TI(:,:,8:end);
    Im = mean(Im, 3);
    srpd = Im ./ Im_TI(:,:,1);
    for ii=1:2  
        v = roi_timeseries(srpd, mask_file, ii, 1);
        p = v.m;
        SR_PD(1, ii) = p(:);      
    end

    % plot LUT

    Gd_LUT = [0:0.01:5];
    
    % with slice profile
    bt = 1.6;
    N = 16;
    [LUT_TI_m, signal_SR, signal_PD] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT', FA_PD, FA_Perf, TD, TR, E1_full, accel_factor, seq_type, seq_type_PD, T1_0_myo(1), T2_0_myo(1), T1_0_myo(1), T2_0_myo(1), r1, r2, bt, N, []);

    % without slice profile
    bt = 1.6;
    N = 1;
    [LUT_TI_m2, signal_SR, signal_PD] = Matlab_gt_perfusion_bloch_simulation(Gd_LUT', FA_PD, FA_Perf, TD, TR, E1_full, accel_factor, seq_type, seq_type_PD, T1_0_myo(1), T2_0_myo(1), T1_0_myo(1), T2_0_myo(1), r1, r2, bt, N, []);
    
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

    params.T1_0 = T1_0_myo(1);
    params.T2_0 = T2_0_myo(1);
    params.T1_0_PD = T1_0_myo(1);
    params.T2_0_PD = T2_0_myo(1);
    
    if(strcmp(seq_type, 'SSFP')==1)
        [LUTslice, Gd] = perfusion_lut_flash_pd_ssfp_sr(params, Gd_LUT);
    else
        [LUTslice, Gd] = perfusion_lut_flash_pd_flash_sr(params, Gd_LUT);
    end
    
    figure
    hold on
    plot(Gd_LUT, LUT_TI_m);
    plot(Gd_LUT, LUT_TI_m2, 'k');
    plot(Gd_LUT, LUTslice(:), 'r');
    hold off
    xlabel('Gd')
    ylabel('SR/PD')
    legend('With slice profile', 'Without slice profile')
    box on
    grid on
    
    sr_over_pd = SR_PD(1, 2);
    perf_Gd_without_slice_profile = interp1(LUTslice, Gd_LUT, sr_over_pd) 
    perf_Gd_with_slice_profile = interp1(LUT_TI_m, Gd_LUT, sr_over_pd)
    perf_Gd_without_slice_profile2 = interp1(LUT_TI_m2, Gd_LUT, sr_over_pd)
end

Gd_perf_slice = perf_Gd_with_slice_profile;
