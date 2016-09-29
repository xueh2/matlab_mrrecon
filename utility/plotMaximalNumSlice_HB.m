
function [aif_duration, perf, TD, Gd, SRcenter, ImagingDuration] = plotMaximalNumSlice_HB(home, sub_dir, TD_in)
% [aif_duration, perf, TD, Gd, SRcenter, ImagingDuration] = plotMaximalNumSlice_HB(home, sub_dir)
% [aif_duration, perf, TD, Gd, SRcenter, ImagingDuration] = plotMaximalNumSlice_HB('E:\gtuser\gt_windows_setup\ut\DualBolus\SSFP', '20150930_A04057_4ml\stress\meas_MID00391_FID82501_BIG__05_STRESS_SSFP_PERF_TPAT3_PF3_4_192x111')

aif_N_runup = 3
aif_seq_type = 'Flash'
aif_Gd_method = 'LUT'

FA_PD = 5
N_runup = 3
seq_type = 'SSFP'
seq_type_PD = 'Flash'
Gd_method = 'LUT'

SRPrep = 33;

ind = find(sub_dir == '\');

v = sub_dir;

filename = v(ind(end)+1:end)

h5Name = fullfile(home, sub_dir, [filename '.h5']);

dset = ismrmrd.Dataset(h5Name);
header = ismrmrd.xml.deserialize(dset.readxml());

FA_Perf = header.sequenceParameters.flipAngle_deg(1)
NumOfPD = header.userParameters.userParameterLong(2).value

if(numel(header.userParameters.userParameterLong)==2)
    NumOfConcatenations = 1
else
    NumOfConcatenations = header.userParameters.userParameterLong(3).value
end

seq_type = header.sequenceParameters.sequence_type
header.acquisitionSystemInformation.systemFieldStrength_T

if(header.acquisitionSystemInformation.systemFieldStrength_T>2)
    r1 = 4.918800;
    r2 = 5.913300;
    if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
        FA_Perf = FA_Perf * ((header.userParameters.userParameterDouble(6).value / header.userParameters.userParameterDouble(5).value) / 0.8436)
    else
        FA_Perf = FA_Perf * ((header.userParameters.userParameterDouble(6).value / header.userParameters.userParameterDouble(5).value) / 0.44433)
    end
else
    r1 = 5.565;
    r2 = 5.78;
end

if(FA_Perf>header.sequenceParameters.flipAngle_deg(1))
    FA_Perf = header.sequenceParameters.flipAngle_deg(1)
end

aif_echo_time_0 = header.sequenceParameters.TE(2)
aif_echo_time_1 = header.sequenceParameters.TE(3)
aif_FA_PD = header.sequenceParameters.flipAngle_deg(2)
aif_FA_Perf = header.sequenceParameters.flipAngle_deg(2)
if(header.sequenceParameters.TI(2)==0)
    aif_TR = 2.12
else
    aif_TR = header.sequenceParameters.echo_spacing(2)
end
aif_E1_full = header.encoding(2).encodingLimits.kspace_encoding_step_1.maximum + 1
aif_accel_factor = header.encoding(2).parallelImaging.accelerationFactor.kspace_encoding_step_1

if(header.sequenceParameters.TI(2)==0)
    aif_TD = 4.7
    aif_TS = aif_TD + ceil( floor(aif_E1_full / aif_accel_factor) / 2.0 ) * aif_TR
else
    aif_TS = header.sequenceParameters.TI(2)
    aif_TD = aif_TS - ceil( floor(aif_E1_full / aif_accel_factor) / 2.0 ) * aif_TR
end

aif_TD = aif_TD + aif_TR

FA_PD = 5
TR = header.sequenceParameters.echo_spacing(1)
E1_full = header.encoding(1).encodedSpace.matrixSize.y
accel_factor = header.encoding(1).parallelImaging.accelerationFactor.kspace_encoding_step_1

if(header.sequenceParameters.TI(2)==0)
    if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
        TD = 17
        TS = TD + TR * ceil( floor(E1_full / accel_factor) / 2.0) - N_runup * TR
    else
        TD = 8.5
        TS = TD + TR * ceil( floor(E1_full / accel_factor) / 2.0);
    end
else
    TS = header.sequenceParameters.TI(1)
    if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
        TD = TS - TR * ceil( floor(E1_full / accel_factor) / 2.0) - N_runup * TR
    else
        TD = TS - TR * ceil( floor(E1_full / accel_factor) / 2.0)
    end
end

perf = SRPrep + TD + ceil( floor(E1_full / accel_factor) * 0.75 ) * TR;
% perf = SRPrep + TS + ceil( floor(E1_full / accel_factor) / 2.0 ) * TR;

if(strcmp(header.sequenceParameters.sequence_type, 'SSFP'))
    perf = perf + N_runup * TR;
end


SLC = header.encoding(1).encodingLimits.slice.maximum+1;

SLCUsed = SLC/NumOfConcatenations;

aif_duration = SRPrep + aif_TS + ceil( floor(aif_E1_full / aif_accel_factor) / 2.0 ) * aif_TR
perf

ImagingDuration = aif_duration + SLCUsed*perf;
60*1e3/ImagingDuration

T1_0_myo = 1300;
T2_0_myo = 50;

params.TR = TR/1e3;
params.TD = TD/1e3;

if(nargin>2)
    params.TD = TD_in/1e3;
    TD = TD_in;
end

params.PD_flip_angle = FA_PD;
params.SR_flip_angle = FA_Perf;
params.offresonance = 0;
params.Npe_full = E1_full;
params.PAT_accel_factor = accel_factor;
params.steadystate_prep = 'linear';
params.N_runup = 3;
params.rf_phase_spoiler_increment = 112;
params.r1 = r1;
params.r2 = r2;
params.T1_0 = T1_0_myo/1e3;
params.T2_0 = T2_0_myo/1e3;
params

Gd = [0:0.01:2];

if(strcmp(seq_type, 'Flash'))
    [LUT_m, Gd, SRcenter, PDcenter] = perfusion_lut_flash_pd_flash_sr(params, Gd');
else
    [LUT_m, Gd, SRcenter, PDcenter] = perfusion_lut_flash_pd_ssfp_sr(params, Gd');
end


