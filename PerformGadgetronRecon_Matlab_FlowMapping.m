
function PerformGadgetronRecon_Matlab_FlowMapping(h5Name, res_dir)
% PerformGadgetronRecon_Matlab_FlowMapping(h5Name, res_dir)

UTDir = getenv('GTPLUS_UT_DIR')

aif_N_runup = 3
aif_seq_type = 'Flash'
aif_Gd_method = 'LUT'

FA_PD = 5
N_runup = 3
seq_type = 'SSFP'
seq_type_PD = 'Flash'
Gd_method = 'LUT'

T1_0_blood_3T = 2000
T2_0_blood_3T = 250
T1_0_myo_3T = 1300
T2_0_myo_3T = 45

T1_0_blood_1p5T = 1700
T2_0_blood_1p5T = 220
T1_0_myo_1p5T = 1100
T2_0_myo_1p5T = 45

%% read in h5 file
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
aif0 = analyze75read('AIF_input_for_moco_0_MAG.hdr');
aif_LUT = analyze75read('aif_cin_LUT_flash_pd_flash_sr.hdr');
aif_PD_echo0 = analyze75read('input_aif_PD_0.hdr');
aif_PD_echo1 = analyze75read('input_aif_PD_1.hdr');
dstAcqTimes = analyze75read('dstAcqTimes_0.hdr');
AIF_AcqTimes = analyze75read('AIF_AcqTimes_0.hdr');
perf_LUT = analyze75read('Perf_T1_Correction_LUT.hdr');

sampleinterval = 0.5;
sigmas = [1.6 4.0 5.3];
sigmaMeasure = 0.1;
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
        mask = analyze75read(fullfile(flow_dir, 'DebugOutput', 'AifLVMask_after_Picking.hdr'));
    catch
        mask = readGTPlusExportImageSeries(flow_dir, 116, 1);
        mask = squeeze(mask);
        mask = mask(:,:,end);
    end
end

figure; imagescn(mask);

numPD = NumOfPD * NumOfConcatenations;
        
for slc=1:SLC

    cd(fullfile(res_dir, 'DebugOutput'))

    Perf_AcqTimes = analyze75read(['Perf_AcqTimes_' num2str(slc-1) '.hdr']);

    % grappa image
    coilmap = readGTPlusExportData(['./slep_2DT_motion_compen_coilMapS_' num2str(slc-1)]);

    RO = size(coilmap, 1)
    E1 = size(coilmap, 2)

    grappaMulCha = readGTPlusExportData(['./slep_2DT_motion_compen_mocoInitial_' num2str(slc-1)]);
    grappaCxIm = SensitivityCoilCombination( grappaMulCha, coilmap);
    grappaCxIm = Zero_Padding_Resize_NoFiltering(grappaCxIm, RO, E1);
    grappa_mag = abs(grappaCxIm);

    % nonlinear image
    nlCxMulCha = readGTPlusExportData(['./slep_2DT_motion_compen_res_final_' num2str(slc-1)]);

    nlCxIm = SensitivityCoilCombination( ifft2c(nlCxMulCha), coilmap);
    nlCxIm = Zero_Padding_Resize_NoFiltering(nlCxIm, RO, E1);
    nl_mag = abs(nlCxIm);

    % -------------------------------
    % perform PSIR recon

    avePDMOCO = mean(nlCxIm(:,:,1:numPD), 3);

    rawFilterRO = generateKSpaceFilter('Hanning', 'Medium', RO, [1 RO], RO/2+1);
    rawFilterE1 = generateKSpaceFilter('Hanning', 'Medium', E1, [1 E1], E1/2+1);

    avePDMOCO_Filtered = performRawDataFilter(fft2c(avePDMOCO), rawFilterRO, rawFilterE1);
    avePDMOCO_Filtered = ifft2c(avePDMOCO_Filtered);
    figure; imagescn(abs(avePDMOCO_Filtered));

    % psir
    backgroundPhs = avePDMOCO_Filtered;
    backgroundPhs = backgroundPhs ./ (abs(backgroundPhs)+eps);

    grappaCxImV = grappaCxIm .* repmat(conj(backgroundPhs), [1 1 size(grappaCxIm, 3)]);
    grappaCxImPSIR = real(grappaCxImV);
    grappaCxImImag = imag(grappaCxImV);

    nlCxImPSIR = nlCxIm .* repmat(conj(backgroundPhs), [1 1 size(nlCxIm, 3)]);
    nlCxImPSIR = real(nlCxImPSIR);

%     grappa = abs(grappaCxImPSIR);
%     nl = abs(nlCxImPSIR);

     grappa = grappaCxImPSIR;
     nl = nlCxImPSIR;

    %% match the scal
    nl2 = nl;

    N = size(grappa, 3);
    for k=1:N
        ak = grappa(:,:,k);
        m1 = norm(ak(:));

        bk = nl(:,:,k);
        m2 = norm(bk(:));

        m1/m2;
        nl2(:,:,k) = nl(:,:,k)* m1/m2;
    end

    nl = nl2;

    figure; imagescn(cat(4, grappa, nl), [], [], [], 3);

    % ----------------------------
    % mag
    nl2 = nl;
    N = size(grappa_mag, 3);
    for k=1:N
        ak = grappa_mag(:,:,k);
        m1 = norm(ak(:));

        bk = nl_mag(:,:,k);
        m2 = norm(bk(:));

        m1/m2;
        nl2(:,:,k) = nl_mag(:,:,k)* m1/m2;
    end

    nl_mag = nl2;

    %% upsample perf

    N = size(grappa, 3)

    perf_PD = squeeze(grappa(:,:,1:numPD,:));
    perf = squeeze(grappa(:,:, numPD+1:N,:));

    % -------------------------------------

    perf_mask = perf_PD(:,:,1);
    ind = find(perf_mask(:)>4);
    perf_mask(:) = 0;
    perf_mask(ind(:)) = 1;
    figure; imagescn(perf_mask);

    RO = size(perf_mask, 1);
    E1 = size(perf_mask, 2);

    ra = 0.20;
    sRO=round(ra*RO);
    eRO = RO - sRO;

    perf_mask(1:sRO, :) = 0;
    perf_mask(eRO:end, :) = 0;

    sE1=round(ra*E1);
    eE1 = E1 - sE1;

    perf_mask(:, 1:sE1) = 0;
    perf_mask(:, eE1:end) = 0;

    figure; imagescn(perf_mask);

    % -------------------------------------

    upsamplingRatio = 2;
    deltaT = 0.5;

    grappa_upsampled = perf;
    nl_upsampled = nl(:,:, numPD+1:end);            

    perf_mag_PD = squeeze(grappa_mag(:,:,1:numPD,:));
    perf_mag = squeeze(grappa_mag(:,:, numPD+1:end,:));

    grappa_mag_upsampled = perf_mag;
    nl_mag_upsampled = nl_mag(:,:, numPD+1:end); 

    %% get the perf SR norm 
    figure; imagescn(cat(4, grappa_upsampled, nl_upsampled), [], [], [], 3);

    N = size(grappa_upsampled, 3);
    PD = perf_PD(:,:,1);

    use_Mask = 1;
    boxFilterSize = 7;
    noisebackground = 1;
    thresRatioForNoise = 7;
    [scc, mask_scc] = Matlab_gt_surface_coil_correction(single(PD), [], use_Mask, boxFilterSize, noisebackground, thresRatioForNoise);

%     PD2 = medfilt2(PD, [3, 3]);
%     SRNorm_Med = grappa_upsampled ./ repmat(PD2+eps, [1 1 N]);

    SRNorm = grappa_upsampled ./ repmat(scc+eps, [1 1 N]);
    SRNorm_NL = nl_upsampled ./ repmat(scc+eps, [1 1 N]);
    figure; imagescn(cat(4, SRNorm, SRNorm_NL), [], [], [], 3);

%             SRNorm_half = grappa_half_upsampled ./ repmat(scc+eps, [1 1 N]);
%             SRNorm_half_NL = nl_half_upsampled ./ repmat(scc+eps, [1 1 N]);
%             figure; imagescn(cat(4, SRNorm_half, SRNorm_half_NL), [], [], [], 3);

    PD_mag = perf_mag_PD(:,:,1);
    [scc, mask_scc] = Matlab_gt_surface_coil_correction(single(PD_mag), [], use_Mask, boxFilterSize, noisebackground, thresRatioForNoise);
    SRNorm_mag = grappa_mag_upsampled ./ repmat(scc+eps, [1 1 N]);
    SRNorm_NL_mag = nl_mag_upsampled ./ repmat(scc+eps, [1 1 N]);

%             SRNorm_half_mag = grappa_half_mag_upsampled ./ repmat(scc+eps, [1 1 N]);
%             SRNorm_half_NL_mag = nl_half_mag_upsampled ./ repmat(scc+eps, [1 1 N]);

    %% convert perf SR norm to Gd
    Gd = 0:0.01:10;
    if(numel(perf_LUT)>3000)
        Gd = 0:0.002:10.002;
    end

    grappa_corrected_noUpSampled = PerformQuantitativePerfusion_GetGdFromLUT(SRNorm, Gd, perf_LUT(1:numel(Gd)));
    nl_corrected_noUpSampled = PerformQuantitativePerfusion_GetGdFromLUT(SRNorm_NL, Gd, perf_LUT(1:numel(Gd)));

    grappa_corrected_mag_noUpSampled = PerformQuantitativePerfusion_GetGdFromLUT(SRNorm_mag, Gd, perf_LUT(1:numel(Gd)));
    nl_corrected_mag_noUpSampled = PerformQuantitativePerfusion_GetGdFromLUT(SRNorm_NL_mag, Gd, perf_LUT(1:numel(Gd)));

    grappa_corrected = PerformQuantitativePerfusion_OneImageSeries_Upsample(grappa_corrected_noUpSampled, dstAcqTimes, Perf_AcqTimes);
    nl_corrected = PerformQuantitativePerfusion_OneImageSeries_Upsample(nl_corrected_noUpSampled, dstAcqTimes, Perf_AcqTimes);

    grappa_corrected_mag = PerformQuantitativePerfusion_OneImageSeries_Upsample(grappa_corrected_mag_noUpSampled, dstAcqTimes, Perf_AcqTimes);
    nl_corrected_mag = PerformQuantitativePerfusion_OneImageSeries_Upsample(nl_corrected_mag_noUpSampled, dstAcqTimes, Perf_AcqTimes);

    %% half sample the gd signal
    grappa_half_corrected = PerformQuantitativePerfusion_OneImageSeries_Upsample(grappa_corrected_noUpSampled(:, :, 1:2:end), dstAcqTimes, Perf_AcqTimes(1:2:end));
    nl_half_corrected = PerformQuantitativePerfusion_OneImageSeries_Upsample(nl_corrected_noUpSampled(:, :, 1:2:end), dstAcqTimes, Perf_AcqTimes(1:2:end));

    grappa_half_corrected_mag = PerformQuantitativePerfusion_OneImageSeries_Upsample(grappa_corrected_mag_noUpSampled(:, :, 1:2:end), dstAcqTimes, Perf_AcqTimes(1:2:end));
    nl_half_corrected_mag = PerformQuantitativePerfusion_OneImageSeries_Upsample(nl_corrected_mag_noUpSampled(:, :, 1:2:end), dstAcqTimes, Perf_AcqTimes(1:2:end));

    figure; imagescn( cat(4, abs(grappa_corrected), abs(nl_corrected), abs(grappa_half_corrected), abs(nl_half_corrected)), [0 1.5], [], [], 3);

    figure; imagescn( cat(4, abs(grappa_corrected), abs(grappa_corrected_mag)), [], [], [], 3);           

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
    DebugFolder = 'E:/gtuser/mrprogs/install/DebugFolder/';
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
        if (reComputed | ~isFileExist(fullfile(home, sub_dir, ['flowmaps_Q_e_' res '_' num2str(slc-1) '.mat'])))

            Fp = [0.1:0.15:6.0001];
            Vp = [0.04:0.01:0.1001];
            PS = [0.35:0.05:4.0001];
            Visf = [0.1:0.075:0.475001];
            numel(Fp)*numel(Vp)*numel(PS)*numel(Visf)

            L     = 0.1;    % cm
            xdelt = L/20;   % cm
            xspan = 0:xdelt:L;

            Gp    = 0;      % ml/g/sec
            Gisf  = 0;      % ml/g/sec
            Dp    = 1e-5;   % cm^2/sec
            Disf  = 1e-6;   % cm^2/sec

            N = numel(aif_cin_Gd_baseline_corrected);

            tic
            [sol_m, C_e_m, C_p_m, Q_e_m] = Matlab_gt_BTEX20_model( double(aif_cin_Gd_baseline_corrected(:)), [0:1:N-1]*deltaT, xspan, Fp, Vp, PS, Visf, Gp, Gisf, Dp, Disf);
            toc

            cd(fullfile(home, sub_dir))
            save(['flowmaps_Q_e_' res '_' num2str(slc-1)], 'Q_e_m', 'Fp', 'Vp', 'PS', 'Visf');
        else
            cd(fullfile(home, sub_dir))
            load(['flowmaps_Q_e_' res '_' num2str(slc-1)]);
        end            
    end

    %% linear recon

    grappa_corrected_with_PSIR = grappa_corrected_mag;
    grappa_corrected_with_PSIR(:,:,1:foot) = grappa_corrected(:,:,1:foot);

    grappa_half_corrected_with_PSIR = grappa_half_corrected_mag;
    grappa_half_corrected_with_PSIR(:,:,1:foot) = grappa_corrected(:,:,1:foot);

    nl_with_PSIR = nl_corrected_mag;
    nl_with_PSIR(:,:,1:foot) = grappa_corrected(:,:,1:foot);

    nl_half_with_PSIR = nl_half_corrected_mag;
    nl_half_with_PSIR(:,:,1:foot) = nl_half_corrected(:,:,1:foot);

    if (reComputed | ~isFileExist(fullfile(home, sub_dir, ['flowmaps_input_' res '_' num2str(slc-1) '.mat'])))
        cd(fullfile(home, sub_dir))
        save(['flowmaps_input_' res '_' num2str(slc-1)], 'grappa_corrected_with_PSIR', 'grappa_half_corrected_with_PSIR', 'nl_with_PSIR', 'nl_half_with_PSIR', 'aif_cin_Gd_baseline_corrected', 'foot', 'peak', 'maxCin');
    end

    if (reComputed | ~isFileExist(fullfile(home, sub_dir, ['flowmaps_Linear_' res '_' suffix '_' num2str(slc-1) '.mat'])))
        tic
        [flowmaps_grappa_PSIR, grappa_interVolumeMap_grappa_PSIR, grappa_MTT_grappa_PSIR, grappa_ecv_grappa_PSIR, Ki_whole_grappa_PSIR, blood_volume_maps_grappa_PSIR, PS_maps_grappa_PSIR, SD_maps_grappa_PSIR] = Matlab_gt_perfusion_flow_mapping(single(aif_cin_Gd_baseline_corrected), single(grappa_corrected_with_PSIR), single(perf_mask), foot, peak-1, maxCin-1, deltaT*1000, data_length_FPWH_ratio, orderBSpline_L1BSpline, numOfInternalControlPoints_L1BSpline, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi, max_iter_TwoCompExp, max_iter_TwoCompFermi, max_iter_BTEX20, max_func_eval_BTEX20, local_search_BTEX20, full_optmization_BTEX20, max_time_shift, fix_shift, hematocrit, two_comp_data_range, compute_BTEX_SD_maps, Q_e_m, Fp, PS, Vp, Visf, DebugFolder);
        toc

        figure; imagescn(cat(3, Ki_whole_grappa_PSIR, flowmaps_grappa_PSIR(:,:,end)), [0 6]); PerfColorMap;
        figure; imagescn(SD_maps_grappa_PSIR(:,:,1), [0 1]); PerfColorMap;
        figure; imagescn(blood_volume_maps_grappa_PSIR(:,:,end), [0 40]); PerfColorMap;

        cd(fullfile(home, sub_dir))
        save(['flowmaps_Linear_' res '_' suffix '_' num2str(slc-1)], ... 
            'flowmaps_grappa_PSIR', 'grappa_interVolumeMap_grappa_PSIR', 'grappa_MTT_grappa_PSIR', 'grappa_ecv_grappa_PSIR', 'Ki_whole_grappa_PSIR', 'blood_volume_maps_grappa_PSIR', 'PS_maps_grappa_PSIR', 'SD_maps_grappa_PSIR');

    end

    if (reComputed | ~isFileExist(fullfile(home, sub_dir, ['flowmaps_Linear_NonLinear_2RR_' res '_'  suffix '_' num2str(slc-1) '.mat'])))
        tic
        [flowmaps_grappa_half_PSIR, grappa_interVolumeMap_grappa_half_PSIR, grappa_MTT_grappa_half_PSIR, grappa_ecv_grappa_half_PSIR, Ki_whole_grappa_half_PSIR, blood_volume_maps_grappa_half_PSIR, PS_maps_grappa_half_PSIR, SD_maps_grappa_half_PSIR] = Matlab_gt_perfusion_flow_mapping(single(aif_cin_Gd_baseline_corrected), single(grappa_half_corrected_with_PSIR), single(perf_mask), foot, peak-1, maxCin-1, deltaT*1000, data_length_FPWH_ratio, orderBSpline_L1BSpline, numOfInternalControlPoints_L1BSpline, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi, max_iter_TwoCompExp, max_iter_TwoCompFermi, max_iter_BTEX20, max_func_eval_BTEX20, local_search_BTEX20, full_optmization_BTEX20, max_time_shift, fix_shift, hematocrit, two_comp_data_range, compute_BTEX_SD_maps, Q_e_m, Fp, PS, Vp, Visf, DebugFolder);
        toc

        cd(fullfile(home, sub_dir))

        save(['flowmaps_Linear_NonLinear_2RR_' res '_' suffix '_' num2str(slc-1)], 'grappa_corrected_with_PSIR', 'grappa_half_corrected_with_PSIR', 'nl_with_PSIR', 'nl_half_with_PSIR', 'aif_cin_Gd_baseline_corrected', 'foot', 'peak', 'maxCin', ... 
            'flowmaps_grappa_half_PSIR',            'grappa_interVolumeMap_grappa_half_PSIR',           'grappa_MTT_grappa_half_PSIR',              'grappa_ecv_grappa_half_PSIR',              'Ki_whole_grappa_half_PSIR',            'blood_volume_maps_grappa_half_PSIR', 'PS_maps_grappa_half_PSIR', 'SD_maps_grappa_half_PSIR');

    end

    if (reComputed | ~isFileExist(fullfile(home, sub_dir, ['flowmaps_Linear_without_R2Star_' res '_'  suffix '_' num2str(slc-1) '.mat'])))
        tic
        [flowmaps_grappa_PSIR_without_R2Star, grappa_interVolumeMap_grappa_PSIR_without_R2Star, grappa_MTT_grappa_PSIR_without_R2Star, grappa_ecv_grappa_PSIR_without_R2Star, Ki_whole_grappa_PSIR_without_R2Star, blood_volume_maps_grappa_PSIR_without_R2Star, PS_maps_grappa_PSIR_without_R2Star, SD_maps_grappa_PSIR_without_R2Star] = Matlab_gt_perfusion_flow_mapping(single(aif_cin_Gd_without_R2Star_baseline_corrected), single(grappa_corrected_with_PSIR), single(perf_mask), foot, peak-1, maxCin-1, deltaT*1000, data_length_FPWH_ratio, orderBSpline_L1BSpline, numOfInternalControlPoints_L1BSpline, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi, max_iter_TwoCompExp, max_iter_TwoCompFermi, max_iter_BTEX20, max_func_eval_BTEX20, local_search_BTEX20, full_optmization_BTEX20, max_time_shift, fix_shift, hematocrit, two_comp_data_range, compute_BTEX_SD_maps, Q_e_m, Fp, PS, Vp, Visf, DebugFolder);
        toc

        cd(fullfile(home, sub_dir))
        save(['flowmaps_Linear_without_R2Star_' res '_' suffix '_' num2str(slc-1)], ... 
            'flowmaps_grappa_PSIR_without_R2Star',  'grappa_interVolumeMap_grappa_PSIR_without_R2Star', 'grappa_MTT_grappa_PSIR_without_R2Star',    'grappa_ecv_grappa_PSIR_without_R2Star',    'Ki_whole_grappa_PSIR_without_R2Star',  'blood_volume_maps_grappa_PSIR_without_R2Star', 'PS_maps_grappa_PSIR_without_R2Star');
    end

    %% nonlinear recon

    if (reComputed | ~isFileExist(fullfile(home, sub_dir, ['flowmaps_NonLinear_' res '_'  suffix '_' num2str(slc-1) '.mat'])))
        tic
        [flowmaps_nl_PSIR, grappa_interVolumeMap_nl_PSIR, grappa_MTT_nl_PSIR, grappa_ecv_nl_PSIR, Ki_whole_nl_PSIR, blood_volume_maps_nl_PSIR, PS_maps_nl_PSIR, SD_maps_nl_PSIR] = Matlab_gt_perfusion_flow_mapping(single(aif_cin_Gd_baseline_corrected), single(nl_with_PSIR), single(perf_mask), foot, peak-1, maxCin-1, deltaT*1000, data_length_FPWH_ratio, orderBSpline_L1BSpline, numOfInternalControlPoints_L1BSpline, max_iter_L1BSpline, lambda_L1BSpline, obj_thres_L1BSpline, grad_thres_L1BSpline, print_iter_L1BSpline, num_of_wavLevels_L1BSpline, with_approx_coeff_L1BSpline, max_iter_Fermi, max_iter_TwoCompExp, max_iter_TwoCompFermi, max_iter_BTEX20, max_func_eval_BTEX20, local_search_BTEX20, full_optmization_BTEX20, max_time_shift, fix_shift, hematocrit, two_comp_data_range, compute_BTEX_SD_maps, Q_e_m, Fp, PS, Vp, Visf, DebugFolder);
        toc

        cd(fullfile(home, sub_dir))
        save(['flowmaps_NonLinear_' res '_' suffix '_' num2str(slc-1)], ... 
            'flowmaps_nl_PSIR', 'grappa_interVolumeMap_nl_PSIR', 'grappa_MTT_nl_PSIR', 'grappa_ecv_nl_PSIR', 'Ki_whole_nl_PSIR', 'blood_volume_maps_nl_PSIR', 'PS_maps_nl_PSIR', 'SD_maps_nl_PSIR');

    end

    closeall

    %% save results

%             save(['flowmaps_Linear_NonLinear_' num2str(slc-1)], 'grappa_corrected_with_PSIR', 'grappa_half_corrected_with_PSIR', 'nl_with_PSIR', 'nl_half_with_PSIR', 'aif_cin_Gd_baseline_corrected', 'foot', 'peak', 'maxCin', ... 
%                 'flowmaps_grappa_PSIR',                 'grappa_interVolumeMap_grappa_PSIR',                'grappa_MTT_grappa_PSIR',                   'grappa_ecv_grappa_PSIR',                   'Ki_whole_grappa_PSIR',                 'blood_volume_maps_grappa_PSIR', 'PS_maps_grappa_PSIR', ...
%                 'flowmaps_grappa_half_PSIR',            'grappa_interVolumeMap_grappa_half_PSIR',           'grappa_MTT_grappa_half_PSIR',              'grappa_ecv_grappa_half_PSIR',              'Ki_whole_grappa_half_PSIR',            'blood_volume_maps_grappa_half_PSIR', 'PS_maps_grappa_half_PSIR', ...
%                 'flowmaps_grappa_PSIR_without_R2Star',  'grappa_interVolumeMap_grappa_PSIR_without_R2Star', 'grappa_MTT_grappa_PSIR_without_R2Star',    'grappa_ecv_grappa_PSIR_without_R2Star',    'Ki_whole_grappa_PSIR_without_R2Star',  'blood_volume_maps_grappa_PSIR_without_R2Star', 'PS_maps_grappa_PSIR_without_R2Star', ...
%                 'flowmaps_nl_PSIR',                     'grappa_interVolumeMap_nl_PSIR',                    'grappa_MTT_nl_PSIR',                       'grappa_ecv_nl_PSIR',                       'Ki_whole_nl_PSIR',                     'blood_volume_maps_nl_PSIR', 'PS_maps_nl_PSIR', ... 
%                 'flowmaps_nl_half_PSIR',                'grappa_interVolumeMap_nl_half_PSIR',               'grappa_MTT_nl_half_PSIR',                  'grappa_ecv_nl_half_PSIR',                  'Ki_whole_nl_half_PSIR',                'blood_volume_maps_nl_half_PSIR', 'PS_maps_nl_half_PSIR');

%             save(['flowmaps_Linear_NonLinear_' num2str(slc-1)], 'grappa_corrected_with_PSIR', 'grappa_half_corrected_with_PSIR', 'nl_with_PSIR', 'nl_half_with_PSIR', 'aif_cin_Gd_baseline_corrected', 'aif_cin_Gd_without_R2Star_baseline_corrected', 'foot', 'peak', 'maxCin', ... 
%                 'flowmaps_grappa_PSIR',                 'grappa_interVolumeMap_grappa_PSIR',                'grappa_MTT_grappa_PSIR',                   'grappa_ecv_grappa_PSIR',                   'Ki_whole_grappa_PSIR',                  ...
%                 'flowmaps_grappa_half_PSIR',            'grappa_interVolumeMap_grappa_half_PSIR',           'grappa_MTT_grappa_half_PSIR',              'grappa_ecv_grappa_half_PSIR',              'Ki_whole_grappa_half_PSIR',             ...
%                 'flowmaps_grappa_PSIR_without_R2Star',  'grappa_interVolumeMap_grappa_PSIR_without_R2Star', 'grappa_MTT_grappa_PSIR_without_R2Star',    'grappa_ecv_grappa_PSIR_without_R2Star',    'Ki_whole_grappa_PSIR_without_R2Star', 'blood_volume_maps_grappa_PSIR_without_R2Star', 'PS_maps_grappa_PSIR_without_R2Star', ...
%                 'flowmaps_nl_PSIR',                     'grappa_interVolumeMap_nl_PSIR',                    'grappa_MTT_nl_PSIR',                       'grappa_ecv_nl_PSIR',                       'Ki_whole_nl_PSIR', ... 
%                 'flowmaps_nl_half_PSIR',                'grappa_interVolumeMap_nl_half_PSIR',               'grappa_MTT_nl_half_PSIR',                  'grappa_ecv_nl_half_PSIR',                  'Ki_whole_nl_half_PSIR');
end
