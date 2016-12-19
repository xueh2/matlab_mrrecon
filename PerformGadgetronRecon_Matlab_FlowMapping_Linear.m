
function PerformGadgetronRecon_Matlab_FlowMapping_Linear(dataDir, h5Name, resDir, roiDir, dataRole, HCT, reComputed, reComputed_OnlyGlobal, reComputed_withoutR2Star)
% PerformGadgetronRecon_Matlab_FlowMapping_Linear(dataDir, h5Name, resDir, roiDir, dataRole, HCT, reComputed)
% PerformGadgetronRecon_Matlab_FlowMapping_Linear('I:\BARTS', 'Perfusion_AIF_TwoEchoes_Interleaved_R2_42110_204543119_204543128_323_20160614-120954', 'I:\ReconResults\BARTS')
% PerformGadgetronRecon_Matlab_FlowMapping_Linear('I:\KAROLINSKA', 'Perfusion_AIF_TwoEchoes_Interleaved_R2_41672_2309137_2309147_949_20160825-093814', 'I:\ReconResults\KAROLINSKA')

if(nargin<4)
    HCT = 0.42;
end

if(nargin<5)
    reComputed = 0;
end

if(nargin<6)
    reComputed_OnlyGlobal = 0;
end

if(nargin<7)
    reComputed_withoutR2Star = 0;
end

UTDir = getenv('GTPLUS_UT_DIR');

aif_N_runup = 3;
aif_seq_type = 'Flash';
aif_Gd_method = 'LUT';

FA_PD = 5;
N_runup = 3;
seq_type = 'SSFP';
seq_type_PD = 'Flash';
Gd_method = 'LUT';

T1_0_blood_3T = 2000;
T2_0_blood_3T = 250;
T1_0_myo_3T = 1300;
T2_0_myo_3T = 45;

T1_0_blood_1p5T = 1700;
T2_0_blood_1p5T = 220;
T1_0_myo_1p5T = 1100;
T2_0_myo_1p5T = 45;

[configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(h5Name);
res_dir = fullfile(resDir, study_dates, h5Name);

res = 'Linear';

hct_str = num2str(HCT);
ind = find(hct_str=='.');
if(~isempty(ind))
    hct_str(ind(:)) = 'p';
end
suffix = [study_dates '_hct' hct_str];

Q_e_name = fullfile(res_dir, ['flowmaps_Q_e_' res '_' suffix '_0.mat']);

slc1_name = fullfile(res_dir, ['flowmaps_Linear_' res '_' suffix '_0.mat']);
slc2_name = fullfile(res_dir, ['flowmaps_Linear_' res '_' suffix '_1.mat']);
slc3_name = fullfile(res_dir, ['flowmaps_Linear_' res '_' suffix '_2.mat']);

processed_linear = 0;
if(~reComputed & isFileExist(Q_e_name) & isFileExist(slc1_name) & isFileExist(slc2_name) & isFileExist(slc3_name))
    disp(['Already processed - ' res_dir ' - ' dataRole ' - hct : ' num2str(HCT)]);
    processed_linear = 1;
end

Q_e_name_without_R2Star = fullfile(res_dir, ['flowmaps_without_R2Star_Q_e_' res '_' suffix '_0.mat']);

slc1_name = fullfile(res_dir, ['flowmaps_Linear_without_R2Star_' res '_' suffix '_0.mat']);
slc2_name = fullfile(res_dir, ['flowmaps_Linear_without_R2Star_' res '_' suffix '_1.mat']);
slc3_name = fullfile(res_dir, ['flowmaps_Linear_without_R2Star_' res '_' suffix '_2.mat']);

processed_linear_withoutR2Star = 0;
if(~reComputed_withoutR2Star & isFileExist(Q_e_name) & isFileExist(slc1_name) & isFileExist(slc2_name) & isFileExist(slc3_name))
    disp(['Already processed without R2* correction - ' res_dir ' - ' dataRole ' - hct : ' num2str(HCT)]);
    processed_linear_withoutR2Star = 1;
end

% if(processed_linear & processed_linear_withoutR2Star)
if(processed_linear)
    return;
end

%% read in h5 file
dset = ismrmrd.Dataset(fullfile(dataDir, [h5Name '.h5']));
header = ismrmrd.xml.deserialize(dset.readxml());

FA_Perf = header.sequenceParameters.flipAngle_deg(1)
NumOfPD = header.userParameters.userParameterLong(2).value

NumOfConcatenations = 1;

for ii=1:numel(header.userParameters.userParameterLong)
    if(strcmp(header.userParameters.userParameterLong(ii).name, 'NumOfConcatenations')==1)
        NumOfConcatenations = header.userParameters.userParameterLong(ii).value;
        break;
    end
    
    if(strcmp(header.userParameters.userParameterLong(ii).name, 'NumOfProtonDensityImages')==1)
        NumOfPD = header.userParameters.userParameterLong(ii).value;
        break;
    end
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

    T1_0_blood = T1_0_blood_3T
    T2_0_blood = T2_0_blood_3T
    T1_0_myo = T1_0_myo_3T
    T2_0_myo = T2_0_myo_3T
else
    r1 = 5.565;
    r2 = 5.78;

    T1_0_blood = T1_0_blood_1p5T
    T2_0_blood = T2_0_blood_1p5T
    T1_0_myo = T1_0_myo_1p5T
    T2_0_myo = T2_0_myo_1p5T
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

SLC = header.encoding(1).encodingLimits.slice.maximum+1;

%% load images

cd(fullfile(res_dir, 'DebugOutput'))

aif_0 = analyze75read('aif_for_TwoEcho_T2StartCorrection_0.hdr');
aif_cin = analyze75read('aif_cin.hdr');
aif_cin_Gd = analyze75read('aif_cin_all_echo0_LUTCorrection.hdr');
aif_cin_Gd_without_R2Star = analyze75read('cin_all_echo0_without_R2Star_LUTCorrection.hdr');
aif_cin_Gd_baseline_corrected = analyze75read('aif_cin_echo0_all_signal_baseline_corrected.hdr');
aif_cin_all_echo0_signal = analyze75read('aif_cin_all_echo0_signal.hdr');
aif_cin_all_echo1_signal = analyze75read('aif_cin_all_echo1_signal.hdr');
aif_cin_all_echo0_signal_after_R2StarCorrection = analyze75read('aif_cin_all_echo0_signal_after_R2StarCorrection.hdr');
aif_cin_all_echo0_OverPD_after_R2StarCorrection = analyze75read('aif_cin_all_echo0_OverPD_after_R2StarCorrection.hdr');
aif_cin_all_R2Star = analyze75read('aif_cin_all_R2Star.hdr');
aif_cin_all_R2Star_SLEP = analyze75read('aif_cin_all_R2Star_SLEP.hdr');
aif_PD = analyze75read('aifPD_for_TwoEcho_T2StartCorrection_0.hdr');
aif_mask = analyze75read('aif_LV_mask_for_TwoEcho_T2StartCorrection_0.hdr');
aif_mask_final = analyze75read('AifLVMask_after_Picking.hdr');
PDMaskIntensity = analyze75read('PDMaskIntensity.hdr');
pdPicked = analyze75read('pdPicked.hdr');
aifPicked = analyze75read('aifPicked.hdr');
aifMaskIntensity = analyze75read('aifMaskIntensity.hdr');
% aif0 = analyze75read('AIF_input_for_moco_0_MAG.hdr');
aif_LUT = analyze75read('aif_cin_LUT_flash_pd_flash_sr.hdr');
% aif_PD_echo0 = analyze75read('input_aif_PD_0.hdr');
% aif_PD_echo1 = analyze75read('input_aif_PD_1.hdr');
dstAcqTimes = analyze75read('dstAcqTimes_0.hdr');
AIF_AcqTimes = analyze75read('AIF_AcqTimes_0.hdr');
perf_LUT = analyze75read('Perf_T1_Correction_LUT.hdr');

sr_0 = analyze75read('input_for_SRNorm_0.hdr');
pd_0 = analyze75read('inputPD_for_SRNorm_0.hdr');
sr_norm_0 = analyze75read('SRNorm_0.hdr');
pd0 = analyze75read('PD_0.hdr');
perf_mask_0 = analyze75read('perf_mask_0.hdr');
gd0 = analyze75read('CASignal_Perf_0.hdr');
Perf_AcqTimes_0 = analyze75read('Perf_AcqTimes_0.hdr');
gd0_resampled = analyze75read('Input_perf_computeFlowMap_0.hdr');

sr_1 = analyze75read('input_for_SRNorm_1.hdr');
pd_1 = analyze75read('inputPD_for_SRNorm_1.hdr');
sr_norm_1 = analyze75read('SRNorm_1.hdr');
pd1 = analyze75read('PD_1.hdr');
perf_mask_1 = analyze75read('perf_mask_1.hdr');
gd1 = analyze75read('CASignal_Perf_1.hdr');
Perf_AcqTimes_1 = analyze75read('Perf_AcqTimes_1.hdr');
gd1_resampled = analyze75read('Input_perf_computeFlowMap_1.hdr');

sr_2 = analyze75read('input_for_SRNorm_2.hdr');
pd_2 = analyze75read('inputPD_for_SRNorm_2.hdr');
sr_norm_2 = analyze75read('SRNorm_2.hdr');
pd2 = analyze75read('PD_2.hdr');
perf_mask_2 = analyze75read('perf_mask_2.hdr');
gd2 = analyze75read('CASignal_Perf_2.hdr');
Perf_AcqTimes_2 = analyze75read('Perf_AcqTimes_2.hdr');
gd2_resampled = analyze75read('Input_perf_computeFlowMap_2.hdr');

% -------------------------------------
RO = size(perf_mask_0, 1);
E1 = size(perf_mask_0, 2);

ra = 0.15;
sRO=round(ra*RO);
eRO = RO - sRO;

sE1=round(ra*E1);
eE1 = E1 - sE1;

pref_prefix = ['perf_mask_' dataRole '_'];

perf_mask_file0 = fullfile(roiDir, [pref_prefix '0.mat']);
if(isFileExist(perf_mask_file0))
    pm = load(perf_mask_file0);
    perf_mask_0 = single(pm.BW);
    perf_mask_0 = flipdim(perf_mask_0, 2);
end

perf_mask_file1 = fullfile(roiDir, [pref_prefix '1.mat']);
if(isFileExist(perf_mask_file1))
    pm = load(perf_mask_file1);
    perf_mask_1 = single(pm.BW);
    perf_mask_1 = flipdim(perf_mask_1, 2);
end

perf_mask_file2 = fullfile(roiDir, [pref_prefix '2.mat']);
if(isFileExist(perf_mask_file2))
    pm = load(perf_mask_file2);
    perf_mask_2 = single(pm.BW);
    perf_mask_2 = flipdim(perf_mask_2, 2);
end


% perf_mask_0(1:sRO, :) = 0;
% perf_mask_0(eRO:end, :) = 0;
% perf_mask_1(1:sRO, :) = 0;
% perf_mask_1(eRO:end, :) = 0;
% perf_mask_2(1:sRO, :) = 0;
% perf_mask_2(eRO:end, :) = 0;
% 
% perf_mask_0(:, 1:sE1) = 0;
% perf_mask_0(:, eE1:end) = 0;
% perf_mask_1(:, 1:sE1) = 0;
% perf_mask_1(:, eE1:end) = 0;
% perf_mask_2(:, 1:sE1) = 0;
% perf_mask_2(:, eE1:end) = 0;

figure; imagescn(cat(3, perf_mask_0, perf_mask_1, perf_mask_2));

% -------------------------------------

gd0_upsampled = PerformQuantitativePerfusion_OneImageSeries_Upsample(gd0, dstAcqTimes, Perf_AcqTimes_0);
gd0_half_upsampled = PerformQuantitativePerfusion_OneImageSeries_Upsample(gd0(:, :, 1:2:end), dstAcqTimes, Perf_AcqTimes_0(1:2:end));

gd1_upsampled = PerformQuantitativePerfusion_OneImageSeries_Upsample(gd1, dstAcqTimes, Perf_AcqTimes_1);
gd1_half_upsampled = PerformQuantitativePerfusion_OneImageSeries_Upsample(gd1(:, :, 1:2:end), dstAcqTimes, Perf_AcqTimes_1(1:2:end));

gd2_upsampled = PerformQuantitativePerfusion_OneImageSeries_Upsample(gd2, dstAcqTimes, Perf_AcqTimes_2);
gd2_half_upsampled = PerformQuantitativePerfusion_OneImageSeries_Upsample(gd2(:, :, 1:2:end), dstAcqTimes, Perf_AcqTimes_2(1:2:end));

figure; imagescn(cat(4, gd0_upsampled, gd0_half_upsampled), [0 2], [], [], 3);

% -------------------------------------

% sampleinterval = 0.5;
% sigmas = [1.6 4.0 5.3];
% sigmaMeasure = 0.1;
% thresGradient = 0.5;

sampleinterval = 0.5;
sigmas = [0.8 2.0 3.2];
sigmaMeasure = 0.2;
thresGradient = 0.5;

plotFlag = 1;

[slope_rest, timeToPeak_rest, peakTime, areaUnderCurve_rest, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(aif_cin_Gd_baseline_corrected(:), sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');

aif_rest2 = aif_cin_Gd_baseline_corrected(round(peakTime):end);
[slope, timeToPeak, peakTime2, areaUnderCurve, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(-aif_rest2(:), sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, 0, 'Feature detection of AIF LV signal');   
valleyTime = round(peakTime2) + peakTime - 0.5;
foot_rest = round(peakTime - timeToPeak_rest);   

v = mean(aif_cin_Gd_without_R2Star(4:foot_rest));
aif_cin_Gd_without_R2Star_baseline_corrected = aif_cin_Gd_without_R2Star - v;

figure
hold on
plot(aif_cin_Gd_baseline_corrected);
plot(aif_cin_Gd_without_R2Star_baseline_corrected, 'r');
hold off

N = size(aif_cin_Gd_baseline_corrected, 1)
NUsed = size(aif_cin, 1)

foot = N-NUsed
peak = round(peakTime)
maxCin = round(valleyTime)

cd(res_dir)

maskFile = fullfile(res_dir, 'DebugOutput', 'aif_mask.mat');

if(isFileExist(maskFile))           

    numAsRep = 0;
    if (NumOfConcatenations>1)
        numAsRep = 1;
    end

    [AIF, acq_time_AIF, physio_time_AIF] = readGTPlusExportImageSeries(flow_dir, 1104, 1, numAsRep);
    [AIF_moco, acq_time_AIF_moco, physio_time_AIF] = readGTPlusExportImageSeries(flow_dir, 101, 1, numAsRep);

    AIF = squeeze(AIF);
    AIF_moco = squeeze(AIF_moco);   
    acq_time_AIF = squeeze(acq_time_AIF);
    acq_time_AIF_moco = squeeze(acq_time_AIF_moco);

    if(numAsRep==1)
        s = size(AIF);
        s(4) = s(4)/2;
        AIF_temp = zeros(s);

        AIF_temp(:,:,1,:) = AIF(:,:,1,1:2:end);
        AIF_temp(:,:,2,:) = AIF(:,:,2,2:2:end);
        AIF = AIF_temp;

        AIF_temp(:,:,1,:) = AIF_moco(:,:,1,1:2:end);
        AIF_temp(:,:,2,:) = AIF_moco(:,:,2,2:2:end);
        AIF_moco = AIF_temp;

        acq_time = zeros(2, s(4));
        acq_time(1,:) = acq_time_AIF_moco(1,1:2:end);
        acq_time(2,:) = acq_time_AIF_moco(2,2:2:end);

        acq_time_AIF_moco = acq_time;
    end

    load (maskFile);
    BW=zeros(size(AIF_moco(:,:,1, 1)));
    BW=roipoly(AIF_moco(:,:,1,1),ROI_info_table(1, 1).ROI_x_coordinates,ROI_info_table(1, 1).ROI_y_coordinates);
    index=find(BW >0);
    mask = BW;
else
    try
        mask = analyze75read(fullfile(res_dir, 'DebugOutput', 'AifLVMask_after_Picking.hdr'));
    catch
        mask = readGTPlusExportImageSeries(res_dir, 116, 1);
        mask = squeeze(mask);
        mask = mask(:,:,end);
    end
end

figure; imagescn(mask);

numPD = NumOfPD * NumOfConcatenations;

if(HCT>0 & HCT<1)
    hct_r = 1/(1-HCT);
    aif_cin_Gd_baseline_corrected = hct_r * aif_cin_Gd_baseline_corrected;
    aif_cin_Gd_without_R2Star_baseline_corrected = hct_r * aif_cin_Gd_without_R2Star_baseline_corrected;
    
    baseV = mean(aif_cin_Gd_baseline_corrected(1:foot));
    aif_cin_Gd_baseline_corrected = aif_cin_Gd_baseline_corrected - baseV;
    
    baseV = mean(aif_cin_Gd_without_R2Star_baseline_corrected(1:foot));
    aif_cin_Gd_without_R2Star_baseline_corrected = aif_cin_Gd_without_R2Star_baseline_corrected - baseV;    
else
    hct_r = 1;
end

for slc=1:SLC

    disp(['========================================================']);
    disp(['Processing slice - ' num2str(slc)]);

    Perf_AcqTimes = analyze75read(fullfile(res_dir, 'DebugOutput', ['Perf_AcqTimes_' num2str(slc-1) '.hdr']));

    %% deconvolution

    data_length_FPWH_ratio = 0;
    orderBSpline_L1BSpline = 3;
    numOfInternalControlPoints_L1BSpline = 7;

    orderBSpline_L1BSpline = 3;
    numOfInternalControlPoints_L1BSpline = 5;        
    max_iter_L1BSpline = 100;
    lambda_L1BSpline = 0.01;
    obj_thres_L1BSpline = 1e-6;
    grad_thres_L1BSpline = 1e-6;
    print_iter_L1BSpline = 0;
    num_of_wavLevels_L1BSpline = 1;
    with_approx_coeff_L1BSpline = 1;
    max_iter_Fermi = 200;
    max_iter_TwoCompExp = 100;
    max_iter_TwoCompFermi = 0;
    max_iter_BTEX20 = 15;
    max_func_eval_BTEX20 = 15;
    hematocrit = 0.45;
    DebugFolder = 'D:/gtuser/mrprogs/install/DebugFolder/';
    deltaT = 0.5;
    data_length_FPWH_ratio = 1;
    max_time_shift = 2.5;
    fix_shift = 0;
    local_search_BTEX20 = 1;
    full_optmization_BTEX20 = 1;
    compute_BTEX_SD_maps = 1;

    two_comp_data_range = 'Whole';

    %% build LUT
    if(slc==1)
        
        Fp = [0.1:0.1:6.0001];
        Vp = [0.01:0.01:0.1001];
        PS = [0.35:0.05:3.001];
        % Visf = [0.1:0.075:0.475001];
        Visf = [0.1:0.075:0.75001];
        disp(['number of Q_e entries - ' num2str(numel(Fp)*numel(Vp)*numel(PS)*numel(Visf))])

        L     = 0.1;    % cm
        xdelt = L/20;   % cm
        xspan = 0:xdelt:L;

        Gp    = 0;      % ml/g/sec
        Gisf  = 0;      % ml/g/sec
        Dp    = 1e-5;   % cm^2/sec
        Disf  = 1e-6;   % cm^2/sec

        N = numel(aif_cin_Gd_baseline_corrected);

            
        if (reComputed | ~isFileExist(Q_e_name))

            if(~isFileExist(Q_e_name))
                disp(['--> compute Q_e_m <--']);
                tic
                [sol_m, C_e_m, C_p_m, Q_e_m] = Matlab_gt_BTEX20_model( double(aif_cin_Gd_baseline_corrected(:)), [0:1:N-1]*deltaT, xspan, Fp, Vp, PS, Visf, Gp, Gisf, Dp, Disf);
                toc

                save(Q_e_name, 'Q_e_m', 'Fp', 'Vp', 'PS', 'Visf');
            else
                load(Q_e_name);
            end
        else
            load(Q_e_name);
        end            
        
        if(0) % if(reComputed_withoutR2Star | ~isFileExist(Q_e_name_without_R2Star))
            disp(['--> compute Q_e_m_without_R2Star <--']);
            tic
            [sol_m_without_R2Star, C_e_m_without_R2Star, C_p_m_without_R2Star, Q_e_m_without_R2Star] = Matlab_gt_BTEX20_model( double(aif_cin_Gd_without_R2Star_baseline_corrected(:)), [0:1:N-1]*deltaT, xspan, Fp, Vp, PS, Visf, Gp, Gisf, Dp, Disf);
            toc

            save(Q_e_name_without_R2Star, 'Q_e_m_without_R2Star', 'Fp', 'Vp', 'PS', 'Visf');
        else
            %load(Q_e_name_without_R2Star);
        end
    end

    %% linear recon

    cd(res_dir)
    
    if(slc==1)
        perf_mask = perf_mask_0;
        grappa_corrected_with_PSIR = gd0_upsampled;
        grappa_half_corrected_with_PSIR = gd0_half_upsampled;
    elseif(slc ==2)
        perf_mask = perf_mask_1;
        grappa_corrected_with_PSIR = gd1_upsampled;
        grappa_half_corrected_with_PSIR = gd1_half_upsampled;
    else
        perf_mask = perf_mask_2;
        grappa_corrected_with_PSIR = gd2_upsampled;
        grappa_half_corrected_with_PSIR = gd2_half_upsampled;
    end
        
    if (reComputed | ~isFileExist(fullfile(res_dir, ['flowmaps_Linear_' res '_' suffix '_' num2str(slc-1) '.mat'])))
        
        disp(['--> compute flow map with HCT = ' hct_str]);
        tic
        [flowmaps_grappa_PSIR, grappa_interVolumeMap_grappa_PSIR, grappa_MTT_grappa_PSIR, grappa_ecv_grappa_PSIR, Ki_whole_grappa_PSIR, blood_volume_maps_grappa_PSIR, PS_maps_grappa_PSIR, SD_maps_grappa_PSIR] = Matlab_gt_perfusion_flow_mapping(single(aif_cin_Gd_baseline_corrected), single(grappa_corrected_with_PSIR), single(perf_mask), foot, peak-1, maxCin-1, deltaT*1000, data_length_FPWH_ratio, orderBSpline_L1BSpline, numOfInternalControlPoints_L1BSpline, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi, max_iter_TwoCompExp, max_iter_TwoCompFermi, max_iter_BTEX20, max_func_eval_BTEX20, local_search_BTEX20, full_optmization_BTEX20, max_time_shift, fix_shift, hematocrit, two_comp_data_range, compute_BTEX_SD_maps, Q_e_m, Fp, PS, Vp, Visf, DebugFolder);
        toc

        flowmaps_grappa_PSIR = flowmaps_grappa_PSIR*hct_r;
        Ki_whole_grappa_PSIR = Ki_whole_grappa_PSIR*hct_r;
        % blood_volume_maps_grappa_PSIR = blood_volume_maps_grappa_PSIR*hct_r;
        SD_maps_grappa_PSIR = SD_maps_grappa_PSIR*hct_r;
        
        figure; imagescn(cat(3, Ki_whole_grappa_PSIR, flowmaps_grappa_PSIR(:,:,end)), [0 8]); PerfColorMap;
        figure; imagescn(SD_maps_grappa_PSIR(:,:,1), [0 1]); PerfColorMap;
        figure; imagescn(blood_volume_maps_grappa_PSIR(:,:,end), [0 40]); PerfColorMap;
        figure; imagescn(grappa_interVolumeMap_grappa_PSIR(:,:,end), [0 80]); PerfColorMap;

        save(['flowmaps_Linear_' res '_' suffix '_' num2str(slc-1)], ... 
            'flowmaps_grappa_PSIR', 'grappa_interVolumeMap_grappa_PSIR', 'grappa_MTT_grappa_PSIR', 'grappa_ecv_grappa_PSIR', 'Ki_whole_grappa_PSIR', 'blood_volume_maps_grappa_PSIR', 'PS_maps_grappa_PSIR', 'SD_maps_grappa_PSIR');
    end
    
%     if (reComputed | ~isFileExist(fullfile(res_dir, ['flowmaps_Linear_HCT_' res '_' suffix '_' num2str(slc-1) '.mat'])))
%         r = 1/(1-HCT);
% 
%         tic
%         [flowmaps_grappa_PSIR_HCT, grappa_interVolumeMap_grappa_PSIR_HCT, grappa_MTT_grappa_PSIR_HCT, grappa_ecv_grappa_PSIR_HCT, Ki_whole_grappa_PSIR_HCT, blood_volume_maps_grappa_PSIR_HCT, PS_maps_grappa_PSIR_HCT, SD_maps_grappa_PSIR_HCT] = Matlab_gt_perfusion_flow_mapping(single(r*aif_cin_Gd_baseline_corrected), single(grappa_corrected_with_PSIR), single(perf_mask), foot, peak-1, maxCin-1, deltaT*1000, data_length_FPWH_ratio, orderBSpline_L1BSpline, numOfInternalControlPoints_L1BSpline, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi, max_iter_TwoCompExp, max_iter_TwoCompFermi, max_iter_BTEX20, max_func_eval_BTEX20, local_search_BTEX20, full_optmization_BTEX20, max_time_shift, fix_shift, hematocrit, two_comp_data_range, compute_BTEX_SD_maps, Q_e_m, Fp, PS, Vp, Visf, DebugFolder);
%         toc
% 
%         flowmaps_grappa_PSIR_HCT = flowmaps_grappa_PSIR_HCT*r;
%         Ki_whole_grappa_PSIR_HCT = Ki_whole_grappa_PSIR_HCT*r;
%         blood_volume_maps_grappa_PSIR_HCT = blood_volume_maps_grappa_PSIR_HCT*r;
%         SD_maps_grappa_PSIR_HCT = SD_maps_grappa_PSIR_HCT*r;
%         
%         figure; imagescn(grappa_interVolumeMap_grappa_PSIR_HCT(:,:,end), [0 80]); PerfColorMap;
%         figure; imagescn(cat(3, Ki_whole_grappa_PSIR_HCT, flowmaps_grappa_PSIR_HCT(:,:,end)), [0 6]); PerfColorMap;
%         
%         figure; imagescn(cat(3, flowmaps_grappa_PSIR(:,:,end), flowmaps_grappa_PSIR_HCT(:,:,end)), [0 6]); PerfColorMap;
%         figure; imagescn(cat(3, grappa_interVolumeMap_grappa_PSIR(:,:,end), grappa_interVolumeMap_grappa_PSIR_HCT(:,:,end)), [0 80]); PerfColorMap;
%         figure; imagescn(cat(3, blood_volume_maps_grappa_PSIR(:,:,end), blood_volume_maps_grappa_PSIR_HCT(:,:,end)), [0 40]); PerfColorMap;
%         
%         save(['flowmaps_Linear_HCT_' res '_' suffix '_' num2str(slc-1)], ... 
%             'flowmaps_grappa_PSIR_HCT', 'grappa_interVolumeMap_grappa_PSIR_HCT', 'grappa_MTT_grappa_PSIR_HCT', 'grappa_ecv_grappa_PSIR_HCT', 'Ki_whole_grappa_PSIR_HCT', 'blood_volume_maps_grappa_PSIR_HCT', 'PS_maps_grappa_PSIR_HCT', 'SD_maps_grappa_PSIR_HCT');
%     end
    
    if(0) % if (reComputed_OnlyGlobal | ~isFileExist(fullfile(res_dir, ['flowmaps_Linear_OnlyGlobalSearch_' res '_' suffix '_' num2str(slc-1) '.mat'])))
        disp(['--> compute flow map with only global search with HCT = ' hct_str]);

        local_search_BTEX20_used = 0;
        
        tic
        [flowmaps_grappa_PSIR_OnlyGlobalSearch, grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, grappa_MTT_grappa_PSIR_OnlyGlobalSearch, grappa_ecv_grappa_PSIR_OnlyGlobalSearch, Ki_whole_grappa_PSIR_OnlyGlobalSearch, blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, PS_maps_grappa_PSIR_OnlyGlobalSearch, SD_maps_grappa_PSIR_OnlyGlobalSearch] = Matlab_gt_perfusion_flow_mapping(single(aif_cin_Gd_baseline_corrected), single(grappa_corrected_with_PSIR), single(perf_mask), foot, peak-1, maxCin-1, deltaT*1000, data_length_FPWH_ratio, orderBSpline_L1BSpline, numOfInternalControlPoints_L1BSpline, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi, max_iter_TwoCompExp, max_iter_TwoCompFermi, max_iter_BTEX20, max_func_eval_BTEX20, local_search_BTEX20_used, full_optmization_BTEX20, max_time_shift, fix_shift, hematocrit, two_comp_data_range, compute_BTEX_SD_maps, Q_e_m, Fp, PS, Vp, Visf, DebugFolder);
        toc

        flowmaps_grappa_PSIR_OnlyGlobalSearch = flowmaps_grappa_PSIR_OnlyGlobalSearch*hct_r;
        Ki_whole_grappa_PSIR_OnlyGlobalSearch = Ki_whole_grappa_PSIR_OnlyGlobalSearch*hct_r;
        % blood_volume_maps_grappa_PSIR_OnlyGlobalSearch = blood_volume_maps_grappa_PSIR_OnlyGlobalSearch*hct_r;
        SD_maps_grappa_PSIR_OnlyGlobalSearch = SD_maps_grappa_PSIR_OnlyGlobalSearch*hct_r;
        
%         figure; imagescn(grappa_interVolumeMap_grappa_PSIR_HCT_OnlyGlobalSearch(:,:,end), [0 80]); PerfColorMap;
%         figure; imagescn(cat(3, Ki_whole_grappa_PSIR_HCT_OnlyGlobalSearch, flowmaps_grappa_PSIR_HCT_OnlyGlobalSearch(:,:,end)), [0 6]); PerfColorMap;
%         
%         figure; imagescn(cat(3, flowmaps_grappa_PSIR(:,:,end), flowmaps_grappa_PSIR_HCT(:,:,end), flowmaps_grappa_PSIR_HCT_OnlyGlobalSearch(:,:,end)), [0 6]); PerfColorMap;
%         figure; imagescn(cat(3, grappa_interVolumeMap_grappa_PSIR(:,:,end), grappa_interVolumeMap_grappa_PSIR_HCT(:,:,end), grappa_interVolumeMap_grappa_PSIR_HCT_OnlyGlobalSearch(:,:,end)), [0 80]); PerfColorMap;
%         figure; imagescn(cat(3, blood_volume_maps_grappa_PSIR(:,:,end), blood_volume_maps_grappa_PSIR_HCT(:,:,end), blood_volume_maps_grappa_PSIR_HCT_OnlyGlobalSearch(:,:,end)), [0 40]); PerfColorMap;

        save(['flowmaps_Linear_OnlyGlobalSearch_' res '_' suffix '_' num2str(slc-1)], ... 
            'flowmaps_grappa_PSIR_OnlyGlobalSearch', 'grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch', 'grappa_MTT_grappa_PSIR_OnlyGlobalSearch', 'grappa_ecv_grappa_PSIR_OnlyGlobalSearch', 'Ki_whole_grappa_PSIR_OnlyGlobalSearch', 'blood_volume_maps_grappa_PSIR_OnlyGlobalSearch', 'PS_maps_grappa_PSIR_OnlyGlobalSearch', 'SD_maps_grappa_PSIR_OnlyGlobalSearch');

    end

%     if (reComputed | ~isFileExist(fullfile(res_dir, ['flowmaps_Linear_NonLinear_2RR_' res '_'  suffix '_' num2str(slc-1) '.mat'])))
%         tic
%         [flowmaps_grappa_half_PSIR, grappa_interVolumeMap_grappa_half_PSIR, grappa_MTT_grappa_half_PSIR, grappa_ecv_grappa_half_PSIR, Ki_whole_grappa_half_PSIR, blood_volume_maps_grappa_half_PSIR, PS_maps_grappa_half_PSIR, SD_maps_grappa_half_PSIR] = Matlab_gt_perfusion_flow_mapping(single(aif_cin_Gd_baseline_corrected), single(grappa_half_corrected_with_PSIR), single(perf_mask), foot, peak-1, maxCin-1, deltaT*1000, data_length_FPWH_ratio, orderBSpline_L1BSpline, numOfInternalControlPoints_L1BSpline, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi, max_iter_TwoCompExp, max_iter_TwoCompFermi, max_iter_BTEX20, max_func_eval_BTEX20, local_search_BTEX20, full_optmization_BTEX20, max_time_shift, fix_shift, hematocrit, two_comp_data_range, compute_BTEX_SD_maps, Q_e_m, Fp, PS, Vp, Visf, DebugFolder);
%         toc
% 
%         flowmaps_grappa_half_PSIR = flowmaps_grappa_half_PSIR*hct_r;
%         Ki_whole_grappa_half_PSIR = Ki_whole_grappa_half_PSIR*hct_r;
%         blood_volume_maps_grappa_half_PSIR = blood_volume_maps_grappa_half_PSIR*hct_r;
%         SD_maps_grappa_half_PSIR = SD_maps_grappa_half_PSIR*hct_r;
%         
%         save(['flowmaps_Linear_NonLinear_2RR_' res '_' suffix '_' num2str(slc-1)], 'grappa_corrected_with_PSIR', 'grappa_half_corrected_with_PSIR', 'nl_with_PSIR', 'nl_half_with_PSIR', 'aif_cin_Gd_baseline_corrected', 'foot', 'peak', 'maxCin', ... 
%             'flowmaps_grappa_half_PSIR',            'grappa_interVolumeMap_grappa_half_PSIR',           'grappa_MTT_grappa_half_PSIR',              'grappa_ecv_grappa_half_PSIR',              'Ki_whole_grappa_half_PSIR',            'blood_volume_maps_grappa_half_PSIR', 'PS_maps_grappa_half_PSIR', 'SD_maps_grappa_half_PSIR');
% 
%     end

    if(0) % if (reComputed_withoutR2Star | ~isFileExist(fullfile(res_dir, ['flowmaps_Linear_without_R2Star_' res '_'  suffix '_' num2str(slc-1) '.mat'])))
        
        disp(['--> compute flow map without R2Star correction and with HCT = ' hct_str]);
        
        tic
        [flowmaps_grappa_PSIR_without_R2Star, grappa_interVolumeMap_grappa_PSIR_without_R2Star, grappa_MTT_grappa_PSIR_without_R2Star, grappa_ecv_grappa_PSIR_without_R2Star, Ki_whole_grappa_PSIR_without_R2Star, blood_volume_maps_grappa_PSIR_without_R2Star, PS_maps_grappa_PSIR_without_R2Star, SD_maps_grappa_PSIR_without_R2Star] = Matlab_gt_perfusion_flow_mapping(single(aif_cin_Gd_without_R2Star_baseline_corrected), single(grappa_corrected_with_PSIR), single(perf_mask), foot, peak-1, maxCin-1, deltaT*1000, data_length_FPWH_ratio, orderBSpline_L1BSpline, numOfInternalControlPoints_L1BSpline, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi, max_iter_TwoCompExp, max_iter_TwoCompFermi, max_iter_BTEX20, max_func_eval_BTEX20, local_search_BTEX20, full_optmization_BTEX20, max_time_shift, fix_shift, hematocrit, two_comp_data_range, compute_BTEX_SD_maps, Q_e_m_without_R2Star, Fp, PS, Vp, Visf, DebugFolder);
        toc

        flowmaps_grappa_PSIR_without_R2Star = flowmaps_grappa_PSIR_without_R2Star*hct_r;
        Ki_whole_grappa_PSIR_without_R2Star = Ki_whole_grappa_PSIR_without_R2Star*hct_r;
        % blood_volume_maps_grappa_PSIR_without_R2Star = blood_volume_maps_grappa_PSIR_without_R2Star*hct_r;
        SD_maps_grappa_PSIR_without_R2Star = SD_maps_grappa_PSIR_without_R2Star*hct_r;
        
        save(['flowmaps_Linear_without_R2Star_' res '_' suffix '_' num2str(slc-1)], ... 
            'flowmaps_grappa_PSIR_without_R2Star',  'grappa_interVolumeMap_grappa_PSIR_without_R2Star', 'grappa_MTT_grappa_PSIR_without_R2Star',    'grappa_ecv_grappa_PSIR_without_R2Star',    'Ki_whole_grappa_PSIR_without_R2Star',  'blood_volume_maps_grappa_PSIR_without_R2Star', 'PS_maps_grappa_PSIR_without_R2Star', 'SD_maps_grappa_PSIR_without_R2Star');
    end

    closeall
end
