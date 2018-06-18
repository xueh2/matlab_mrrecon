
function [ROITable, sf, rf, sf_i, rf_i, res_table] = PerformGadgetronRecon_SavedIsmrmrd_Matlab_ROIValues(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, fixed_HCT, report_column, prefix, processing_always, processing_snr_always, processing_Gd_always, prefix_SNR)
% [ROITable, sf, rf, sf_i, rf_i, res_table] = PerformGadgetronRecon_SavedIsmrmrd_Matlab_ROIValues(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, fixed_HCT, report_column, prefix, processing_always)
% if isempty(fixed_HCT)==1, measured HCT is used; if(fixed_HCT>=0 & fixed_HCT<=1), this value is used

if(nargin<11)
    processing_always = 1;
end

if(nargin<12)
    processing_snr_always = 1;
end

if(nargin<13)
    processing_Gd_always = 1;
end

if(nargin<14)
    prefix_SNR = prefix;
end

ROITable = PerfTable(1, :);
ROITable = [ROITable 'stress flow' 'rest flow' 'stress ischemia flow' 'rest ischemia flow' 'stress aif peak Gd with baseline correction' 'rest aif peak Gd with baseline correction' 'hematocrit' ...
    'stress E' 'rest E' 'stress PS' 'rest PS' ... 
    'stress Visf' 'rest Visf' 'stress Vp' 'rest Vp' ... 
    'stress Ki_MF' 'rest Ki_MF' 'stress Ki_Fermi' 'rest Ki_Fermi' ... 
    'stress Ki_TwoCompExp' 'rest Ki_TwoCompExp' 'stress Ki_BTEX' 'rest Ki_BTEX' ... 
    'stress E_i' 'rest E_i' 'stress PS_i' 'rest PS_i' ... 
    'stress Visf_i' 'rest Visf_i' 'stress Vp_i' 'rest Vp_i' ... 
    'stress Ki_MF_i' 'rest Ki_MF_i' 'stress Ki_Fermi_i' 'rest Ki_Fermi_i' ... 
    'stress Ki_TwoCompExp_i' 'rest Ki_TwoCompExp_i' 'stress Ki_BTEX_i' 'rest Ki_BTEX_i' ];

nV = numel(PerfTable(1, :));

scanInd = [];
patientID = [];
scanDate = [];
scanTime = [];
age = [];
gender = [];
height = [];
weight = [];
bsa = [];
diabetes = [];
HTN = [];
Chol = [];
smoker = [];
EF = [];
stressHB = [];
restHB = [];
hematocrit = [];
LVMass = [];
LVEDV = [];

sf = [];
rf = [];
sf_i = [];
rf_i = [];

sE = [];
rE = [];
sE_i = [];
rE_i = [];

sPS = [];
rPS = [];
sPS_i = [];
rPS_i = [];

sVisf = [];
rVisf = [];
sVisf_i = [];
rVisf_i = [];

sVb = [];
rVb = [];
sVb_i = [];
rVb_i = [];

sVp = [];
rVp = [];
sVp_i = [];
rVp_i = [];

sKi_MF = [];
rKi_MF = [];
sKi_MF_i = [];
rKi_MF_i = [];

sKi_Fermi = [];
rKi_Fermi = [];
sKi_Fermi_i = [];
rKi_Fermi_i = [];

sKi_TwoCompExp = [];
rKi_TwoCompExp = [];
sKi_TwoCompExp_i = [];
rKi_TwoCompExp_i = [];

sKi_BTEX = [];
rKi_BTEX = [];
sKi_BTEX_i = [];
rKi_BTEX_i = [];

sSD = [];
rSD = [];
sSD_i = [];
rSD_i = [];

sSNR = [];
rSNR = [];

sGd = [];
rGd = [];

sT2S = [];
rT2S = [];

sAifPeakGd = [];
rAifPeakGd = [];

sSplenic = [];
rSplenic = [];

sDelay = [];
rDelay = [];

sFermiDelay = [];
rFermiDelay = [];

pre_T1_blood = [];
pre_T1_myo = [];

post_T1_blood = [];
post_T1_myo = [];

num_column = size(PerfTable, 2);
num = size(PerfTable, 1)-1;
for n=1:num
       
    stressCase = PerfTable{n+1, stress_column};
    restCase = PerfTable{n+1, rest_column};
    if(ischemia_column<=num_column)
        Category = PerfTable{n+1, ischemia_column};
    else
        Category = 1;
    end
    
    disp([num2str(n) ' out of ' num2str(num) ' - ' PerfTable{n+1, report_column} ' - Processing : '  ' - ' PerfTable{n+1, rest_column}]); 
    disp(['==================================================================']);  
       
    [configName, scannerID, pID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(stressCase);
    [configName, scannerID, pID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(restCase);

    % figDir = fullfile(resDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_Figure']) 
    figDir = fullfile(resDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_' study_time_stress '_' study_time_rest '_Figure']);
    
    roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_ROI'])
    
    stressDir = fullfile(resDir, study_dates, stressCase)
    restDir = fullfile(resDir, study_dates, restCase)
    
    if(~isempty(fixed_HCT))
        if(fixed_HCT>=0 & fixed_HCT<=1)
            HCT = fixed_HCT;
        end
    else
        HCT = PerfTable{n+1, hct_column};
        if(isnan(HCT))
            HCT = 0;
        else
            if(HCT>1)
                HCT = HCT/100;
            end
        end
    end    
    
    disp([num2str(n) ' out of ' num2str(num) ' - HCT : ' num2str(HCT)]); 
   
    hct_str = num2str(HCT);
    ind = find(hct_str=='.');
    if(~isempty(ind))
        hct_str(ind(:)) = 'p';
    end
    suffix = [study_dates '_hct' hct_str];
    
    if(isFileExist(fullfile(roiDir, 'myo_stress1.mat')) | isFileExist(fullfile(roiDir, 'myo_stress2.mat')) | isFileExist(fullfile(roiDir, 'myo_stress3.mat')))
        s1_roi = 'myo_stress1.mat';
        s2_roi = 'myo_stress2.mat';
        s3_roi = 'myo_stress3.mat';
        
        r1_roi = 'myo_rest1.mat';
        r2_roi = 'myo_rest2.mat';
        r3_roi = 'myo_rest3.mat';
    else
        s1_roi = 's1.mat';
        s2_roi = 's2.mat';
        s3_roi = 's3.mat';
        
        r1_roi = 'r1.mat';
        r2_roi = 'r2.mat';
        r3_roi = 'r3.mat';
    end
    
    stress_has_roi = (isFileExist(fullfile(roiDir, s1_roi))) | (isFileExist(fullfile(roiDir, s2_roi))) | (isFileExist(fullfile(roiDir, s3_roi)));
    rest_has_roi = (isFileExist(fullfile(roiDir, r1_roi))) | (isFileExist(fullfile(roiDir, r2_roi))) | (isFileExist(fullfile(roiDir, r3_roi)));
    
    if(stress_has_roi & rest_has_roi)
        
        scanInd = [scanInd; n];
        patientID = [patientID; {pID}];
        scanDate = [scanDate; study_dates];
        scanTime = [scanTime; study_time_stress];
        hematocrit = [hematocrit; HCT];

        age = [age; PerfTable{n+1, 12}];
        gender = [gender; PerfTable{n+1, 13}];

        stressHB = [stressHB; PerfTable{n+1, 8}];
        restHB = [restHB; PerfTable{n+1, 11}];

        height = [height; get_valid_value(PerfTable{n+1, 14})];
        weight = [weight; get_valid_value(PerfTable{n+1, 15})];

        diabetes = [diabetes; get_valid_value(PerfTable{n+1, 17})];
        HTN = [HTN; get_valid_value(PerfTable{n+1, 18})];
        Chol = [Chol; get_valid_value(PerfTable{n+1, 19})];
        smoker = [smoker; get_valid_value(PerfTable{n+1, 20})];
        EF = [EF; get_valid_value(PerfTable{n+1, 27})];
        
        vv = get_valid_value(PerfTable{n+1, 35});
        if(isempty(get_valid_value(PerfTable{n+1, 35})))
            LVMass = [LVMass; -1];
            LVEDV = [LVEDV; -1];
        else            
            LVMass = [LVMass; get_valid_value(PerfTable{n+1, 35})];
            LVEDV = [LVEDV; get_valid_value(PerfTable{n+1, 36})];
        end
        
        disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]);       
        
        two_ROI = 0;

        % PerfResult_file = fullfile(roiDir, [prefix '_PerfResult_' suffix '.mat']);
        PerfResult_file = fullfile(figDir, [prefix '_PerfResult_' suffix '.mat']);
        
        SNRResult_file = fullfile(roiDir, [prefix_SNR '_SNR_PerfResult_' scannerID '_' pID '_' studyID '_' study_dates '.mat']);
        Gd_Result_file = fullfile(roiDir, [prefix_SNR '_Gd_PerfResult_' scannerID '_' pID '_' studyID '_' study_dates '.mat']);
        
%         SNRResult_file = fullfile(figDir, [prefix_SNR '_SNR_PerfResult_' scannerID '_' pID '_' studyID '_' study_dates '.mat']);
%         Gd_Result_file = fullfile(figDir, [prefix_SNR '_Gd_PerfResult_' scannerID '_' pID '_' studyID '_' study_dates '.mat']);
        
        has_result_file = 0;
        if(~processing_always & isFileExist(PerfResult_file) & isFileExist(SNRResult_file) & isFileExist(Gd_Result_file))
            tt = load(PerfResult_file);
            v = tt.PerfResult;
            has_result_file = 1;
            
            res_stress = tt.res_stress;
            res_rest = tt.res_rest;
            two_ROI = tt.two_ROI;
            
            sE = [sE; res_stress.E];
            rE = [rE; res_rest.E];

            sPS = [sPS; res_stress.PS];
            rPS = [rPS; res_rest.PS];

            sVisf = [sVisf; res_stress.Visf];
            rVisf = [rVisf; res_rest.Visf];

            sVp = [sVp; res_stress.Vp];
            rVp = [rVp; res_rest.Vp];

            sVb = [sVb; res_stress.Vb];
            rVb = [rVb; res_rest.Vb];
            
            sKi_MF = [sKi_MF; res_stress.Ki_MF];
            rKi_MF = [rKi_MF; res_rest.Ki_MF];

            sKi_Fermi = [sKi_Fermi; res_stress.Ki_Fermi];
            rKi_Fermi = [rKi_Fermi; res_rest.Ki_Fermi];

            sKi_TwoCompExp = [sKi_TwoCompExp; res_stress.Ki_TwoCompExp];
            rKi_TwoCompExp = [rKi_TwoCompExp; res_rest.Ki_TwoCompExp];

            sKi_BTEX = [sKi_BTEX; res_stress.Ki_BTEX];
            rKi_BTEX = [rKi_BTEX; res_rest.Ki_BTEX];

            sSD = [sSD; res_stress.SD];
            rSD = [rSD; res_rest.SD];
            
            sSplenic = [sSplenic; mean(tt.sSplenic)];
            rSplenic = [rSplenic; mean(tt.rSplenic)];
            
            pre_T1_blood = [pre_T1_blood; tt.pre_T1_blood(end)];
            pre_T1_myo = [pre_T1_myo; tt.pre_T1_myo(end)];
            
            post_T1_blood = [post_T1_blood; tt.post_T1_blood(end)];
            post_T1_myo = [post_T1_myo; tt.post_T1_myo(end)];
            
            sSNR = [sSNR; mean(tt.s_SNR(4:end, :), 1)];
            rSNR = [rSNR; mean(tt.r_SNR(4:end, :), 1)];
            
            sDelay = [sDelay; res_stress.delay];
            rDelay = [rDelay; res_rest.delay];
            
            sFermiDelay = [sFermiDelay; res_stress.fermi_delay];
            rFermiDelay = [rFermiDelay; res_rest.fermi_delay];
            
            if(two_ROI)
                sE_i = [sE_i; res_stress.E_i];
                rE_i = [rE_i; res_rest.E_i];

                sPS_i = [sPS_i; res_stress.PS_i];
                rPS_i = [rPS_i; res_rest.PS_i];

                sVisf_i = [sVisf_i; res_stress.Visf_i];
                rVisf_i = [rVisf_i; res_rest.Visf_i];

                sVp_i = [sVp_i; res_stress.Vp_i];
                rVp_i = [rVp_i; res_rest.Vp_i];
                
                sVb_i = [sVb_i; res_stress.Vb_i];
                rVb_i = [rVb_i; res_rest.Vb_i];

                sKi_MF_i = [sKi_MF_i; res_stress.Ki_MF_i];
                rKi_MF_i = [rKi_MF_i; res_rest.Ki_MF_i];

                sKi_Fermi_i = [sKi_Fermi_i; res_stress.Ki_Fermi_i];
                rKi_Fermi_i = [rKi_Fermi_i; res_rest.Ki_Fermi_i];

                sKi_TwoCompExp_i = [sKi_TwoCompExp_i; res_stress.Ki_TwoCompExp_i];
                rKi_TwoCompExp_i = [rKi_TwoCompExp_i; res_rest.Ki_TwoCompExp_i];

                sKi_BTEX_i = [sKi_BTEX_i; res_stress.Ki_BTEX_i];
                rKi_BTEX_i = [rKi_BTEX_i; res_rest.Ki_BTEX_i];

                sSD_i = [sSD_i; res_stress.SD_i];
                rSD_i = [rSD_i; res_rest.SD_i];
            else
                sE_i = [sE_i; -1 -1 -1];
                rE_i = [rE_i; -1 -1 -1];

                sPS_i = [sPS_i; -1 -1 -1];
                rPS_i = [rPS_i; -1 -1 -1];

                sVisf_i = [sVisf_i; -1 -1 -1];
                rVisf_i = [rVisf_i; -1 -1 -1];

                sVp_i = [sVp_i; -1 -1 -1];
                rVp_i = [rVp_i; -1 -1 -1];
                
                sVb_i = [sVb_i; -1 -1 -1];
                rVb_i = [rVb_i; -1 -1 -1];

                sKi_MF_i = [sKi_MF_i; -1 -1 -1];
                rKi_MF_i = [rKi_MF_i; -1 -1 -1];

                sKi_Fermi_i = [sKi_Fermi_i; -1 -1 -1];
                rKi_Fermi_i = [rKi_Fermi_i; -1 -1 -1];

                sKi_TwoCompExp_i = [sKi_TwoCompExp_i; -1 -1 -1];
                rKi_TwoCompExp_i = [rKi_TwoCompExp_i; -1 -1 -1];

                sKi_BTEX_i = [sKi_BTEX_i; -1 -1 -1];
                rKi_BTEX_i = [rKi_BTEX_i; -1 -1 -1];

                sSD_i = [sSD_i; -1 -1 -1];
                rSD_i = [rSD_i; -1 -1 -1];
            end
            
            snr_v = load(SNRResult_file);
                
            s_SNR = snr_v.s_SNR;
            r_SNR = snr_v.r_SNR;

%             sSNR = [sSNR; max(s_SNR(4:end, :))];
%             rSNR = [rSNR; max(r_SNR(4:end, :))];
            
            gd_v = load(Gd_Result_file);
                
            s_Gd = gd_v.s_Gd;
            r_Gd = gd_v.r_Gd;

            sGd = [sGd; max(s_Gd(1:end, :))];
            rGd = [rGd; max(r_Gd(1:end, :))];
            
            sT2S = [sT2S; max(gd_v.sR2Star(1:end, :))];
            rT2S = [rT2S; max(gd_v.rR2Star(1:end, :))];
            
            if(isfield(gd_v, 'stress_aif_Gd_peak'))
                sAifPeakGd = [sAifPeakGd; gd_v.stress_aif_Gd_peak];
                rAifPeakGd = [rAifPeakGd; gd_v.rest_aif_Gd_peak];
            else
                sAifPeakGd = [sAifPeakGd; -1];
                rAifPeakGd = [rAifPeakGd; -1];
            end
        else        
            cd(roiDir)
            v = PerfTable(n+1, :);

            nV = numel(v);

            if(isFileExist(fullfile(roiDir, s1_roi)))
                s1 = load(fullfile(roiDir, s1_roi));
            else
                s1 = [];
            end
            
            if(isFileExist(fullfile(roiDir, s2_roi)))
                s2 = load(fullfile(roiDir, s2_roi));
            else
                s2 = [];
            end
            
            
            if(isFileExist(fullfile(roiDir, s3_roi)))
                s3 = load(fullfile(roiDir, s3_roi));
            else
                s3 = [];
            end

            if(isFileExist(fullfile(roiDir, r1_roi)))
                r1 = load(fullfile(roiDir, r1_roi));
            else
                r1 = [];
            end 
            
            if(isFileExist(fullfile(roiDir, r2_roi)))
                r2 = load(fullfile(roiDir, r2_roi));
            else
                r2 = [];
            end
            
            if(isFileExist(fullfile(roiDir, r3_roi)))
                r3 = load(fullfile(roiDir, r3_roi));
            else
                r3 = [];
            end

%             if(~isFileExist(fullfile(roiDir, 'rest.mat')))
%                 copyfile(fullfile(figDir, 'rest.mat'), roiDir);
%             end
%             
%             if(~isFileExist(fullfile(roiDir, 'stress.mat')))
%                 copyfile(fullfile(figDir, 'stress.mat'), roiDir);
%             end
            
            rest = load(fullfile(figDir, 'rest.mat'));
            stress = load(fullfile(figDir, 'stress.mat'));
            
            disp(['--> loading mat results : <--']);
            sMat1 = fullfile(resDir, study_dates, stressCase, [prefix '_' suffix '_0.mat'])
            sMat2 = fullfile(resDir, study_dates, stressCase, [prefix '_' suffix '_1.mat'])
            sMat3 = fullfile(resDir, study_dates, stressCase, [prefix '_' suffix '_2.mat'])
            disp(['--> ============================================ <--']);
            
            stressMat1 = load(sMat1);
            stressMat2 = load(sMat2);
            stressMat3 = load(sMat3);
            
            try
                rMat1 = fullfile(resDir, study_dates, restCase, [prefix '_' suffix '_0.mat'])
                restMat1 = load(rMat1);
                restMat2 = load(fullfile(resDir, study_dates, restCase, [prefix '_' suffix '_1.mat']));
                restMat3 = load(fullfile(resDir, study_dates, restCase, [prefix '_' suffix '_2.mat']));
            catch
                restMat1 = stressMat1;
                restMat2 = stressMat2;
                restMat3 = stressMat3;
            end
            
%             figure; imagescn(cat(3, stressMat1.flowmaps_grappa_PSIR(:,:,end), stressMat1.Ki_whole_grappa_PSIR), [0 8]); PerfColorMap;
%             figure; imagescn(cat(3, stressMat2.flowmaps_grappa_PSIR(:,:,end), stressMat2.Ki_whole_grappa_PSIR), [0 8]); PerfColorMap;
%             figure; imagescn(cat(3, stressMat3.flowmaps_grappa_PSIR(:,:,end), stressMat3.Ki_whole_grappa_PSIR), [0 8]); PerfColorMap;
% 
% %             pause;
%             closeall

            try
                btex_res_0 = analyze75read(fullfile(stressDir, 'slc0', 'DebugFolder', 'BTEX_res_0'));
                btex_res_1 = analyze75read(fullfile(stressDir, 'slc1', 'DebugFolder', 'BTEX_res_0'));
                btex_res_2 = analyze75read(fullfile(stressDir, 'slc2', 'DebugFolder', 'BTEX_res_0'));

                fermi_res_0 = analyze75read(fullfile(stressDir, 'slc0', 'DebugFolder', 'Fermi_res_0'));
                fermi_res_1 = analyze75read(fullfile(stressDir, 'slc1', 'DebugFolder', 'Fermi_res_0'));
                fermi_res_2 = analyze75read(fullfile(stressDir, 'slc2', 'DebugFolder', 'Fermi_res_0'));

    %             btex_res_0 = flipdim(btex_res_0, 2);
    %             btex_res_1 = flipdim(btex_res_1, 2);
    %             btex_res_2 = flipdim(btex_res_2, 2);
    %             
    %             fermi_res_0 = flipdim(fermi_res_0, 2);
    %             fermi_res_1 = flipdim(fermi_res_1, 2);
    %             fermi_res_2 = flipdim(fermi_res_2, 2);

                stressMat1.btex_delay = btex_res_0(:,:,end); 
                stressMat2.btex_delay = btex_res_1(:,:,end); 
                stressMat3.btex_delay = btex_res_2(:,:,end); 

                stressMat1.fermi_delay = fermi_res_0(:,:,end); 
                stressMat2.fermi_delay = fermi_res_1(:,:,end); 
                stressMat3.fermi_delay = fermi_res_2(:,:,end); 

                btex_res_0 = analyze75read(fullfile(restDir, 'slc0', 'DebugFolder', 'BTEX_res_0'));
                btex_res_1 = analyze75read(fullfile(restDir, 'slc1', 'DebugFolder', 'BTEX_res_0'));
                btex_res_2 = analyze75read(fullfile(restDir, 'slc2', 'DebugFolder', 'BTEX_res_0'));

                fermi_res_0 = analyze75read(fullfile(restDir, 'slc0', 'DebugFolder', 'Fermi_res_0'));
                fermi_res_1 = analyze75read(fullfile(restDir, 'slc1', 'DebugFolder', 'Fermi_res_0'));
                fermi_res_2 = analyze75read(fullfile(restDir, 'slc2', 'DebugFolder', 'Fermi_res_0'));

                restMat1.btex_delay = btex_res_0(:,:,end); 
                restMat2.btex_delay = btex_res_1(:,:,end); 
                restMat3.btex_delay = btex_res_2(:,:,end); 

                restMat1.fermi_delay = fermi_res_0(:,:,end); 
                restMat2.fermi_delay = fermi_res_1(:,:,end); 
                restMat3.fermi_delay = fermi_res_2(:,:,end); 
            catch
                stressMat1.btex_delay = [];
                stressMat2.btex_delay = [];
                stressMat3.btex_delay = [];
                
                stressMat1.fermi_delay = [];
                stressMat2.fermi_delay = [];
                stressMat3.fermi_delay = [];
                
                restMat1.btex_delay = [];
                restMat2.btex_delay = [];
                restMat3.btex_delay = [];
                
                restMat1.fermi_delay = [];
                restMat2.fermi_delay = [];
                restMat3.fermi_delay = [];
            end
            
            res_stress = PerformGadgetronRecon_Matlab_ROIValues_OneCase(s1, s2, s3, stressMat1, stressMat2, stressMat3, HCT);
            res_rest = PerformGadgetronRecon_Matlab_ROIValues_OneCase(r1, r2, r3, restMat1, restMat2, restMat3, HCT);                       
            
            v{nV+1} = res_stress.flow;
            disp(['stress flow : ' num2str(v{nV+1})]);
            disp(['stress Ki_Fermi : ' num2str(res_stress.Ki_Fermi)]);
            
            v{nV+2} = res_rest.flow;
            disp(['rest flow : ' num2str(v{nV+2})]);

            has_rest_second_roi = 0;
            if( (~isempty(s1) & numel(s1.ROI_info_table)==2) | (~isempty(s2) & numel(s2.ROI_info_table)==2) | (~isempty(s3) & numel(s3.ROI_info_table)==2))
                
                two_ROI = 1;

                v{nV+3} = res_stress.flow_i;
                disp(['stress flow, ischemia : ' num2str(v{nV+3})]);
                
                if(isfield(res_rest, 'flow_i'))
                    v{nV+4} = res_rest.flow_i;
                    has_rest_second_roi = 1;
                else
                    v{nV+4} = [-1 -1 -1];
                end
            else
                v{nV+3} = [-1 -1 -1];
                v{nV+4} = [-1 -1 -1];
            end

            v{nV+5} = max(stress.aif_stress_baseline_corrected);
            v{nV+6} = max(rest.aif_rest_baseline_corrected);
            v{nV+7} = HCT;
            
            ind_s = find(res_stress.flow>0);
            % if(~isempty(find(res_stress.flow(ind_s)>6)) | ~isempty(find(res_stress.flow(ind_s)<3)))
            % if(~isempty(find(res_stress.flow(ind_s)>6)) | ~isempty(find(res_stress.flow(ind_s)<1.5)))
%             if( ~isempty(find(res_rest.flow(ind_s)>1.5)) | ~isempty(find(res_stress.flow(ind_s)>6)) )
%                 if(isfield(stressMat1, 'flowmaps_grappa_PSIR'))
%                     figure; imagescn(flipdim(stressMat1.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, s1_roi)); PerfColorMap;
%                     figure; imagescn(flipdim(stressMat2.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, s2_roi)); PerfColorMap;
%                     figure; imagescn(flipdim(stressMat3.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, s3_roi)); PerfColorMap;
% 
%                     figure; imagescn(flipdim(restMat1.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, r1_roi)); PerfColorMap;
%                     figure; imagescn(flipdim(restMat2.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, r2_roi)); PerfColorMap;
%                     figure; imagescn(flipdim(restMat3.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, r3_roi)); PerfColorMap;
% 
%                     onlyReview = 1;
%                     [h_flow_stress, h_flow_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir,  stressCase, restCase, [0 6], onlyReview, resDir);        
% 
%                     pause;
%                     closeall
%                 end
%             end
            
            if(isfield(stressMat1, 'flowmaps_grappa_PSIR'))
                    if(~isempty(s1)) figure; imagescn(flipdim(stressMat1.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, s1_roi)); PerfColorMap; end
                    if(~isempty(s2)) figure; imagescn(flipdim(stressMat2.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, s2_roi)); PerfColorMap; end
                    if(~isempty(s3)) figure; imagescn(flipdim(stressMat3.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, s3_roi)); PerfColorMap; end

                    if(~isempty(r1)) figure; imagescn(flipdim(restMat1.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, r1_roi)); PerfColorMap; end
                    if(~isempty(r2)) figure; imagescn(flipdim(restMat2.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, r2_roi)); PerfColorMap; end
                    if(~isempty(r3)) figure; imagescn(flipdim(restMat3.flowmaps_grappa_PSIR(:,:,end), 2), [0 8], [], [], [], fullfile(roiDir, r3_roi)); PerfColorMap; end

%                     if(pause_cases) 
%                         user_in = input('accept cases y or n :');
%                         if(user_in=='n')
%                             onlyReview = 1;
%                             [h_flow_stress, h_flow_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir,  stressCase, restCase, [0 6], onlyReview, resDir);      
%                             pause;
%                         end
%                     end
                    closeall
            end
                
            % ---------------------------------------------
            % load splenic
            sp1 = fullfile(roiDir, 's1_splenic.mat');
            sp2 = fullfile(roiDir, 's2_splenic.mat');
            sp3 = fullfile(roiDir, 's3_splenic.mat');
            
            rp1 = fullfile(roiDir, 'r1_splenic.mat');
            rp2 = fullfile(roiDir, 'r2_splenic.mat');
            rp3 = fullfile(roiDir, 'r3_splenic.mat');
            
            stress_splenic = 0;
            rest_splenic = 0;
            
            if(isFileExist(sp1))
                vv = load(sp1);
                stress_splenic = vv.ROI_info_table.ROI_mean;
            elseif(isFileExist(sp2))
                vv = load(sp2);
                stress_splenic = vv.ROI_info_table.ROI_mean;
            elseif(isFileExist(sp3))
                vv = load(sp3);
                stress_splenic = vv.ROI_info_table.ROI_mean;
            end
            
            if(isFileExist(rp1))
                vv = load(rp1);
                rest_splenic = vv.ROI_info_table.ROI_mean;
            elseif(isFileExist(rp2))
                vv = load(rp2);
                rest_splenic = vv.ROI_info_table.ROI_mean;
            elseif(isFileExist(rp3))
                vv = load(rp3);
                rest_splenic = vv.ROI_info_table.ROI_mean;
            end
            
            sSplenic = [sSplenic; stress_splenic];
            rSplenic = [rSplenic; rest_splenic];
            
            % ---------------------------------------------
            % load T1 map
            t1_pre = fullfile(roiDir, 'T1_pre.mat');
            t1_post = fullfile(roiDir, 'T1_post.mat');
            
            pre_t1_blood = [0 0 0];
            pre_t1_myo = [0 0 0];
            post_t1_blood = [0 0 0];
            post_t1_myo = [0 0 0];
            
            if(isFileExist(t1_pre))
                t1_v = load(t1_pre);
                pre_t1_blood = [t1_v.ROI_info_table(1,1).ROI_mean t1_v.ROI_info_table(3,2).ROI_mean t1_v.ROI_info_table(5,3).ROI_mean];
                pre_t1_myo = [t1_v.ROI_info_table(2,1).ROI_mean t1_v.ROI_info_table(4,2).ROI_mean t1_v.ROI_info_table(6,3).ROI_mean];
            end
            
            if(isFileExist(t1_post))
                t1_v = load(t1_post);
                post_t1_blood = [t1_v.ROI_info_table(1,1).ROI_mean t1_v.ROI_info_table(3,2).ROI_mean t1_v.ROI_info_table(5,3).ROI_mean];
                post_t1_myo = [t1_v.ROI_info_table(2,1).ROI_mean t1_v.ROI_info_table(4,2).ROI_mean t1_v.ROI_info_table(6,3).ROI_mean];
            end
            
            pre_T1_blood = [pre_T1_blood; pre_t1_blood];
            pre_T1_myo = [pre_T1_myo; pre_t1_myo];
            
            post_T1_blood = [post_T1_blood; post_t1_blood];
            post_T1_myo = [post_T1_myo; post_t1_myo];

            % ---------------------------------------------
            % SNR
            
            %if(processing_snr_always | ~isFileExist(SNRResult_file))
            if(0)
                try
                    sdata1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'input_spatiotemporal_filter__row0_MAG.hdr'));
                    sdata2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'input_spatiotemporal_filter__row1_MAG.hdr'));
                    sdata3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'input_spatiotemporal_filter__row2_MAG.hdr'));

                    cd(stressDir)
                    s_gfactor = readGTPlusExportImageSeries_Squeeze(300);
                    % s_gfactor = flipdim(s_gfactor, 2);

%                     moco_perf = readGTPlusExportImageSeries_Squeeze(104);
%                     moco_perf = flipdim(moco_perf, 2);
%                     
%                     sdata1 = squeeze(moco_perf(:,:,1,:));
%                     sdata2 = squeeze(moco_perf(:,:,2,:));
%                     sdata3 = squeeze(moco_perf(:,:,3,:));

                    rdata1 = analyze75read(fullfile(restDir, 'DebugOutput', 'input_spatiotemporal_filter__row0_MAG.hdr'));
                    rdata2 = analyze75read(fullfile(restDir, 'DebugOutput', 'input_spatiotemporal_filter__row1_MAG.hdr'));
                    rdata3 = analyze75read(fullfile(restDir, 'DebugOutput', 'input_spatiotemporal_filter__row2_MAG.hdr'));

                    cd(restDir)
                    r_gfactor = readGTPlusExportImageSeries_Squeeze(300);
%                     r_gfactor = flipdim(r_gfactor, 2);

%                     moco_perf = readGTPlusExportImageSeries_Squeeze(104);
%                     moco_perf = flipdim(moco_perf, 2);
%                     
%                     rdata1 = squeeze(moco_perf(:,:,1,:));
%                     rdata2 = squeeze(moco_perf(:,:,2,:));
%                     rdata3 = squeeze(moco_perf(:,:,3,:));
                    
%                     sdata1 = flipdim(sdata1, 2);
%                     sdata2 = flipdim(sdata2, 2);
%                     sdata3 = flipdim(sdata3, 2);
%                     rdata1 = flipdim(rdata1, 2);
%                     rdata2 = flipdim(rdata2, 2);
%                     rdata3 = flipdim(rdata3, 2);

                    num_stress = size(sdata1, 3);
                    snr_s1 = 25 * sdata1 ./ squeeze(s_gfactor(:,:,1,1:num_stress));
                    snr_s2 = 25 * sdata2 ./ squeeze(s_gfactor(:,:,2,1:num_stress));
                    snr_s3 = 25 * sdata3 ./ squeeze(s_gfactor(:,:,3,1:num_stress));

                    num_rest = size(rdata1, 3);
                    snr_r1 = 25 * rdata1 ./ squeeze(r_gfactor(:,:,1,1:num_rest));
                    snr_r2 = 25 * rdata2 ./ squeeze(r_gfactor(:,:,2,1:num_rest));
                    snr_r3 = 25 * rdata3 ./ squeeze(r_gfactor(:,:,3,1:num_rest));

                    cd(roiDir)
                    if(~isempty(s1))
                        BW1=zeros(size(sdata1(:,:,1)));
                        BW1=roipoly(sdata1(:,:,1), s1(1).ROI_info_table(1).ROI_x_coordinates, s1(1).ROI_info_table(1).ROI_y_coordinates);
                        index1=find(BW1 >0);
                    end

                    if(~isempty(s2))
                        BW2=zeros(size(sdata2(:,:,1)));
                        BW2=roipoly(sdata2(:,:,1), s2(1).ROI_info_table(1).ROI_x_coordinates, s2(1).ROI_info_table(1).ROI_y_coordinates);
                        index2=find(BW2 >0);
                    end

                    if(~isempty(s3))
                        BW3=zeros(size(sdata3(:,:,1)));
                        BW3=roipoly(sdata3(:,:,1), s3(1).ROI_info_table(1).ROI_x_coordinates, s3(1).ROI_info_table(1).ROI_y_coordinates);
                        index3=find(BW3 >0);
                    end

                    if(~isempty(r1))
                        rBW1=zeros(size(rdata1(:,:,1)));
                        rBW1=roipoly(rdata1(:,:,1), r1(1).ROI_info_table(1).ROI_x_coordinates, r1(1).ROI_info_table(1).ROI_y_coordinates);
                        rindex1=find(rBW1 >0);
                    end

                    if(~isempty(r2))
                        rBW2=zeros(size(rdata2(:,:,1)));
                        rBW2=roipoly(rdata2(:,:,1), r2(1).ROI_info_table(1).ROI_x_coordinates, r2(1).ROI_info_table(1).ROI_y_coordinates);
                        rindex2=find(rBW2 >0);
                    end

                    if(~isempty(r3))
                        rBW3=zeros(size(rdata3(:,:,1)));
                        rBW3=roipoly(rdata3(:,:,1), r3(1).ROI_info_table(1).ROI_x_coordinates, r3(1).ROI_info_table(1).ROI_y_coordinates);
                        rindex3=find(rBW3 >0);
                    end

                    % gfactor scaled by 100
                    % data is scaled by 4
                    nRep = size(sdata1, 3);

                    s_SNR = zeros(nRep, 3);
                    r_SNR = zeros(nRep, 3);

                    for rr=1:nRep
                        sd1 = snr_s1(:,:,rr);
                        sd2 = snr_s2(:,:,rr);
                        sd3 = snr_s3(:,:,rr);

                        rd1 = snr_r1(:,:,rr);
                        rd2 = snr_r2(:,:,rr);
                        rd3 = snr_r3(:,:,rr);

                        v1 = 0;
                        if(~isempty(s1))
                            v1 = mean(sd1(index1));
                        end

                        v2 = 0;
                        if(~isempty(s2))
                            v2 = mean(sd2(index2));
                        end

                        v3 = 0;
                        if(~isempty(s3))
                            v3 = mean(sd3(index3));
                        end

                        s_SNR(rr, :) = [v1 v2 v3];

                        v1 = 0;
                        if(~isempty(r1))
                            v1 = mean(rd1(rindex1));
                        end

                        v2 = 0;
                        if(~isempty(r2))
                            v2 = mean(rd2(rindex2));
                        end

                        v3 = 0;
                        if(~isempty(r3))
                            v3 = mean(rd3(rindex3));
                        end
                        r_SNR(rr, :) = [v1 v2 v3];                
                    end

                    sSNR = [sSNR; max(s_SNR(4:end, :))];
                    rSNR = [rSNR; max(r_SNR(4:end, :))];

                    s_gfactor_map = s_gfactor(:,:,:,1) /100;
                    r_gfactor_map = r_gfactor(:,:,:,1) /100;
                    save(SNRResult_file, 's_SNR', 'r_SNR', 'snr_s1', 'snr_s2', 'snr_s3', 'snr_r1', 'snr_r2', 'snr_r3', 's_gfactor_map', 'r_gfactor_map');               
                catch
                    s_SNR = 0;
                    r_SNR = 0;
                    sSNR = [sSNR; -1 -1 -1];
                    rSNR = [rSNR; -1 -1 -1];
                end
            else
                try
                snr_v = load(SNRResult_file);
                
                s_SNR = snr_v.s_SNR;
                r_SNR = snr_v.r_SNR;
                
                sSNR = [sSNR; max(s_SNR(4:end, :))];
                rSNR = [rSNR; max(r_SNR(4:end, :))];
                catch
                    sSNR = [sSNR; -1 -1 -1];
                    rSNR = [rSNR; -1 -1 -1];
                end
            end

            % ---------------------------------------------
            % Gd concentration
            if(processing_Gd_always | ~isFileExist(Gd_Result_file))
                try
                    sdata1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_PSIR_0.hdr'));
                    sdata2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_PSIR_1.hdr'));
                    sdata3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_PSIR_2.hdr'));

                    rdata1 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_PSIR_0.hdr'));
                    rdata2 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_PSIR_1.hdr'));
                    rdata3 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_PSIR_2.hdr'));

                    sdata1 = flipdim(sdata1, 2);
                    sdata2 = flipdim(sdata2, 2);
                    sdata3 = flipdim(sdata3, 2);
                    rdata1 = flipdim(rdata1, 2);
                    rdata2 = flipdim(rdata2, 2);
                    rdata3 = flipdim(rdata3, 2);

                    [HeartRate_rest, aif_cin_Gd_rest, aif_cin_Gd_rest_without_R2Star, aif_cin_all_echo0_signal_rest, aif_cin_all_echo0_signal_after_R2StarCorrection_rest, footTime_rest, peakTime_rest, valleyTime_rest, R2Star_rest, KiMap, flowMap, EMap, PSMap, VisfMap, VpMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, restCase);
                    [HeartRate_stress, aif_cin_Gd_stress, aif_cin_Gd_stress_without_R2Star, aif_cin_all_echo0_signal_stress, aif_cin_all_echo0_signal_after_R2StarCorrection_stress, footTime_stress, peakTime_stress, valleyTime_stress, R2Star_stress,KiMap, flowMap, EMap, PSMap, VisfMap, VpMap] = PerformGadgetronRecon_Statistics_PerfusionCase_OneScan(resDir, stressCase);

                    sR2Star = analyze75read(fullfile(stressDir, 'DebugOutput', 'aif_cin_all_R2Star_SLEP.hdr'));
                    rR2Star = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_cin_all_R2Star_SLEP.hdr'));

                    rest_T2S = 1.0 ./ (rR2Star+eps);
                    stress_T2S = 1.0 ./ (sR2Star+eps);

                    rest_aif_T2S_peak = min( [rest_T2S(ceil(peakTime_rest)) rest_T2S(ceil(peakTime_rest+0.5)) rest_T2S(ceil(peakTime_rest-0.5))] );
                    stress_aif_T2S_peak = min( [stress_T2S(ceil(peakTime_stress)) stress_T2S(ceil(peakTime_stress+0.5)) stress_T2S(ceil(peakTime_stress-0.5))] );

                    rest_aif_Gd_peak = max( aif_cin_Gd_rest );
                    stress_aif_Gd_peak = max( aif_cin_Gd_stress );

                    sAifPeakGd = [sAifPeakGd; stress_aif_Gd_peak];
                    rAifPeakGd = [rAifPeakGd; rest_aif_Gd_peak];

                    cd(roiDir)
                    if(~isempty(s1))
                        BW1=zeros(size(sdata1(:,:,1)));
                        BW1=roipoly(sdata1(:,:,1), s1(1).ROI_info_table(1).ROI_x_coordinates, s1(1).ROI_info_table(1).ROI_y_coordinates);
                        index1=find(BW1 >0);
                    end

                    if(~isempty(s2))
                        BW2=zeros(size(sdata2(:,:,1)));
                        BW2=roipoly(sdata2(:,:,1), s2(1).ROI_info_table(1).ROI_x_coordinates, s2(1).ROI_info_table(1).ROI_y_coordinates);
                        index2=find(BW2 >0);
                    end

                    if(~isempty(s3))
                        BW3=zeros(size(sdata3(:,:,1)));
                        BW3=roipoly(sdata3(:,:,1), s3(1).ROI_info_table(1).ROI_x_coordinates, s3(1).ROI_info_table(1).ROI_y_coordinates);
                        index3=find(BW3 >0);
                    end

                    if(~isempty(r1))
                        rBW1=zeros(size(rdata1(:,:,1)));
                        rBW1=roipoly(rdata1(:,:,1), r1(1).ROI_info_table(1).ROI_x_coordinates, r1(1).ROI_info_table(1).ROI_y_coordinates);
                        rindex1=find(rBW1 >0);
                    end

                    if(~isempty(r2))
                        rBW2=zeros(size(rdata2(:,:,1)));
                        rBW2=roipoly(rdata2(:,:,1), r2(1).ROI_info_table(1).ROI_x_coordinates, r2(1).ROI_info_table(1).ROI_y_coordinates);
                        rindex2=find(rBW2 >0);
                    end

                    if(~isempty(r3))
                        rBW3=zeros(size(rdata3(:,:,1)));
                        rBW3=roipoly(rdata3(:,:,1), r3(1).ROI_info_table(1).ROI_x_coordinates, r3(1).ROI_info_table(1).ROI_y_coordinates);
                        rindex3=find(rBW3 >0);
                    end

                    nRep = size(sdata1, 3);

                    s_Gd = zeros(nRep, 3);
                    r_Gd = zeros(nRep, 3);

                    for rr=1:nRep
                        sd1 = sdata1(:,:,rr);
                        sd2 = sdata2(:,:,rr);
                        sd3 = sdata3(:,:,rr);

                        rd1 = rdata1(:,:,rr);
                        rd2 = rdata2(:,:,rr);
                        rd3 = rdata3(:,:,rr);

                        v1 = 0;
                        if(~isempty(s1))
                            v1 = mean(sd1(index1));
                        end

                        v2 = 0;
                        if(~isempty(s2))
                            v2 = mean(sd2(index2));
                        end

                        v3 = 0;
                        if(~isempty(s3))
                            v3 = mean(sd3(index3));
                        end

                        s_Gd(rr, :) = [v1 v2 v3];

                        v1 = 0;
                        if(~isempty(r1))
                            v1 = mean(rd1(rindex1));
                        end

                        v2 = 0;
                        if(~isempty(r2))
                            v2 = mean(rd2(rindex2));
                        end

                        v3 = 0;
                        if(~isempty(r3))
                            v3 = mean(rd3(rindex3));
                        end
                        r_Gd(rr, :) = [v1 v2 v3];                
                    end
                   
                    sGd = [sGd; max(s_Gd(1:end, :))];
                    rGd = [rGd; max(r_Gd(1:end, :))];

%                     vs = mean( [s_Gd(ceil(peakTime_rest)) s_Gd(ceil(peakTime_rest+0.5)) s_Gd(ceil(peakTime_rest-0.5))] );
%                     sGd = [sGd; vs];
%                     
%                     vr = mean( [r_Gd(ceil(peakTime_rest)) r_Gd(ceil(peakTime_rest+0.5)) r_Gd(ceil(peakTime_rest-0.5))] );
%                     rGd = [rGd; vr];

%                     sR2S = [sR2S; max(sR2Star(1:end, :))];
%                     rR2S = [rR2S; max(rR2Star(1:end, :))];

                    sT2S = [sT2S; stress_aif_T2S_peak];
                    rT2S = [rT2S; rest_aif_T2S_peak];

                    save(Gd_Result_file, 's_Gd', 'r_Gd', 'sR2Star', 'rR2Star', 'peakTime_rest', 'peakTime_stress', 'aif_cin_Gd_rest', 'aif_cin_Gd_stress', 'stress_aif_Gd_peak', 'rest_aif_Gd_peak');  
                catch
                    sGd = [sGd; -1 -1 -1];
                    rGd = [rGd; -1 -1 -1];

                    sT2S = [sT2S; -1 -1 -1];
                    rT2S = [rT2S; -1 -1 -1];
                end
            else
                try
                    gd_v = load(Gd_Result_file);

                    s_Gd = gd_v.s_Gd;
                    r_Gd = gd_v.r_Gd;

                    sGd = [sGd; max(s_Gd(1:end, :))];
                    rGd = [rGd; max(r_Gd(1:end, :))];

                    sT2S = [sT2S; max(gd_v.sR2Star(1:end, :))];
                    rT2S = [rT2S; max(gd_v.rR2Star(1:end, :))];

                    if(isfield(gd_v, 'stress_aif_Gd_peak'))
                        sAifPeakGd = [sAifPeakGd; gd_v.stress_aif_Gd_peak];
                        rAifPeakGd = [rAifPeakGd; gd_v.rest_aif_Gd_peak];
                    else
                        sAifPeakGd = [sAifPeakGd; 0];
                        rAifPeakGd = [rAifPeakGd; 0];
                    end
                catch
                    sGd = [sGd; -1 -1 -1];
                    rGd = [rGd; -1 -1 -1];

                    sT2S = [sT2S; -1 -1 -1];
                    rT2S = [rT2S; -1 -1 -1];
                end
            end
                
            % ---------------------------------------------            
            ind = 8;
            
            % E
            v{nV+ind} = res_stress.E; ind = ind+1;
            v{nV+ind} = res_rest.E; ind = ind+1;
            
            % PS
            v{nV+ind} = res_stress.PS; ind = ind+1;
            v{nV+ind} = res_rest.PS; ind = ind+1;
            
            % Visf
            v{nV+ind} = res_stress.Visf; ind = ind+1;
            v{nV+ind} = res_rest.Visf; ind = ind+1;
            
            % Vp
            v{nV+ind} = res_stress.Vp; ind = ind+1;
            v{nV+ind} = res_rest.Vp; ind = ind+1;
            
            % Ki_MF
            v{nV+ind} = res_stress.Ki_MF; ind = ind+1;
            v{nV+ind} = res_rest.Ki_MF; ind = ind+1;
                        
            % Ki_Fermi
            v{nV+ind} = res_stress.Ki_Fermi; ind = ind+1;
            v{nV+ind} = res_rest.Ki_Fermi; ind = ind+1;
            
            % Ki_TwoCompExp
            v{nV+ind} = res_stress.Ki_TwoCompExp; ind = ind+1;
            v{nV+ind} = res_rest.Ki_TwoCompExp; ind = ind+1;
            
            % Ki_BTEX
            v{nV+ind} = res_stress.Ki_BTEX; ind = ind+1;
            v{nV+ind} = res_rest.Ki_BTEX; ind = ind+1;
            
            sE = [sE; res_stress.E];
            rE = [rE; res_rest.E];

            sPS = [sPS; res_stress.PS];
            rPS = [rPS; res_rest.PS];

            sVisf = [sVisf; res_stress.Visf];
            rVisf = [rVisf; res_rest.Visf];

            sVp = [sVp; res_stress.Vp];
            rVp = [rVp; res_rest.Vp];
            
            sVb = [sVb; res_stress.Vb];
            rVb = [rVb; res_rest.Vb];

            sKi_MF = [sKi_MF; res_stress.Ki_MF];
            rKi_MF = [rKi_MF; res_rest.Ki_MF];

            sKi_Fermi = [sKi_Fermi; res_stress.Ki_Fermi];
            rKi_Fermi = [rKi_Fermi; res_rest.Ki_Fermi];

            sKi_TwoCompExp = [sKi_TwoCompExp; res_stress.Ki_TwoCompExp];
            rKi_TwoCompExp = [rKi_TwoCompExp; res_rest.Ki_TwoCompExp];

            sKi_BTEX = [sKi_BTEX; res_stress.Ki_BTEX];
            rKi_BTEX = [rKi_BTEX; res_rest.Ki_BTEX];

            sSD = [sSD; res_stress.SD];
            rSD = [rSD; res_rest.SD];
            
            sDelay = [sDelay; res_stress.delay];
            rDelay = [rDelay; res_rest.delay];
            
            sFermiDelay = [sFermiDelay; res_stress.fermi_delay];
            rFermiDelay = [rFermiDelay; res_rest.fermi_delay];
            
            if(two_ROI)
                
                if(~has_rest_second_roi)
                    res_rest.flow_i = [-1 -1 -1];
                    res_rest.E_i = [-1 -1 -1];
                    res_rest.PS_i = [-1 -1 -1];
                    res_rest.Visf_i = [-1 -1 -1];
                    res_rest.Vp_i = [-1 -1 -1];
                    res_rest.Vb_i = [-1 -1 -1];
                    res_rest.Ki_MF_i = [-1 -1 -1];
                    res_rest.Ki_Fermi_i = [-1 -1 -1];
                    res_rest.Ki_TwoCompExp_i = [-1 -1 -1];
                    res_rest.Ki_BTEX_i = [-1 -1 -1];
                    res_rest.SD_i = [-1 -1 -1];
                end
                
                % E
                v{nV+ind} = res_stress.E_i; ind = ind+1;
                v{nV+ind} = res_rest.E_i; ind = ind+1;

                % PS
                v{nV+ind} = res_stress.PS_i; ind = ind+1;
                v{nV+ind} = res_rest.PS_i; ind = ind+1;

                % Visf
                v{nV+ind} = res_stress.Visf_i; ind = ind+1;
                v{nV+ind} = res_rest.Visf_i; ind = ind+1;

                % Vp
                v{nV+ind} = res_stress.Vp_i; ind = ind+1;
                v{nV+ind} = res_rest.Vp_i; ind = ind+1;

                % Ki_MF
                v{nV+ind} = res_stress.Ki_MF_i; ind = ind+1;
                v{nV+ind} = res_rest.Ki_MF_i; ind = ind+1;

                % Ki_Fermi
                v{nV+ind} = res_stress.Ki_Fermi_i; ind = ind+1;
                v{nV+ind} = res_rest.Ki_Fermi_i; ind = ind+1;

                % Ki_TwoCompExp
                v{nV+ind} = res_stress.Ki_TwoCompExp_i; ind = ind+1;
                v{nV+ind} = res_rest.Ki_TwoCompExp_i; ind = ind+1;

                % Ki_BTEX
                v{nV+ind} = res_stress.Ki_BTEX_i; ind = ind+1;
                v{nV+ind} = res_rest.Ki_BTEX_i; ind = ind+1;
                
                sE_i = [sE_i; res_stress.E_i];
                rE_i = [rE_i; res_rest.E_i];

                sPS_i = [sPS_i; res_stress.PS_i];
                rPS_i = [rPS_i; res_rest.PS_i];

                sVisf_i = [sVisf_i; res_stress.Visf_i];
                rVisf_i = [rVisf_i; res_rest.Visf_i];

                sVp_i = [sVp_i; res_stress.Vp_i];
                rVp_i = [rVp_i; res_rest.Vp_i];
                
                sVb_i = [sVb_i; res_stress.Vb_i];
                rVb_i = [rVb_i; res_rest.Vb_i];

                sKi_MF_i = [sKi_MF_i; res_stress.Ki_MF_i];
                rKi_MF_i = [rKi_MF_i; res_rest.Ki_MF_i];

                sKi_Fermi_i = [sKi_Fermi_i; res_stress.Ki_Fermi_i];
                rKi_Fermi_i = [rKi_Fermi_i; res_rest.Ki_Fermi_i];

                sKi_TwoCompExp_i = [sKi_TwoCompExp_i; res_stress.Ki_TwoCompExp_i];
                rKi_TwoCompExp_i = [rKi_TwoCompExp_i; res_rest.Ki_TwoCompExp_i];

                sKi_BTEX_i = [sKi_BTEX_i; res_stress.Ki_BTEX_i];
                rKi_BTEX_i = [rKi_BTEX_i; res_rest.Ki_BTEX_i];

                sSD_i = [sSD_i; res_stress.SD_i];
                rSD_i = [rSD_i; res_rest.SD_i];
            else
                % E
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % PS
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Visf
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Vp
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Ki_MF
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Ki_Fermi
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Ki_TwoCompExp
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;

                % Ki_BTEX
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                v{nV+ind} = [-1 -1 -1]; ind = ind+1;
                
                sE_i = [sE_i; -1 -1 -1];
                rE_i = [rE_i; -1 -1 -1];

                sPS_i = [sPS_i; -1 -1 -1];
                rPS_i = [rPS_i; -1 -1 -1];

                sVisf_i = [sVisf_i; -1 -1 -1];
                rVisf_i = [rVisf_i; -1 -1 -1];

                sVp_i = [sVp_i; -1 -1 -1];
                rVp_i = [rVp_i; -1 -1 -1];
                
                sVb_i = [sVb_i; -1 -1 -1];
                rVb_i = [rVb_i; -1 -1 -1];

                sKi_MF_i = [sKi_MF_i; -1 -1 -1];
                rKi_MF_i = [rKi_MF_i; -1 -1 -1];

                sKi_Fermi_i = [sKi_Fermi_i; -1 -1 -1];
                rKi_Fermi_i = [rKi_Fermi_i; -1 -1 -1];

                sKi_TwoCompExp_i = [sKi_TwoCompExp_i; -1 -1 -1];
                rKi_TwoCompExp_i = [rKi_TwoCompExp_i; -1 -1 -1];

                sKi_BTEX_i = [sKi_BTEX_i; -1 -1 -1];
                rKi_BTEX_i = [rKi_BTEX_i; -1 -1 -1];

                sSD_i = [sSD_i; -1 -1 -1];
                rSD_i = [rSD_i; -1 -1 -1];
            end
            
            PerfResult = v;
            save(PerfResult_file, 'PerfResult', 'res_stress', 'res_rest', 'two_ROI', 'sSplenic', 'rSplenic', 'pre_T1_blood', 'pre_T1_myo', 'post_T1_blood', 'post_T1_myo');           
        end
    
        ROITable = [ROITable; v];
                
        rf = [rf; ROITable{end, nV+2}];         
        
        if(two_ROI)
            sf_i = [sf_i; ROITable{end, nV+3}];
            rf_i = [rf_i; ROITable{end, nV+4}];
        else
            sf_i = [sf_i; -1 -1 -1];
            rf_i = [rf_i; -1 -1 -1];
        end
        sf = [sf; ROITable{end, nV+1}];       
    else
        disp(['Cannot find the ROI ... ']);
        winopen(figDir);
        mkdir(roiDir);
        winopen(roiDir);
        pause;
    end
end

sf_mean = [];
rf_mean = [];

sE_mean = [];
rE_mean = [];

sVp_mean = [];
rVp_mean = [];

sVb_mean = [];
rVb_mean = [];

sVisf_mean = [];
rVisf_mean = [];

sPS_mean = [];
rPS_mean = [];

sKi_MF_mean = [];
rKi_MF_mean = [];
sKi_Fermi_mean = [];
rKi_Fermi_mean = [];
sKi_TwoCompExp_mean = [];
rKi_TwoCompExp_mean = [];
sKi_BTEX_mean = [];
rKi_BTEX_mean = [];

sDelay_mean = [];
rDelay_mean = [];

sFermiDelay_mean = [];
rFermiDelay_mean = [];

for n=1:size(sf,1)
    sf_mean = [sf_mean; get_entry_mean(sf(n,:))];
    rf_mean = [rf_mean; get_entry_mean(rf(n,:))];

    sE_mean = [sE_mean; get_entry_mean(sE(n,:))];
    rE_mean = [rE_mean; get_entry_mean(rE(n,:))];
    
    sPS_mean = [sPS_mean; get_entry_mean(sPS(n,:))];
    rPS_mean = [rPS_mean; get_entry_mean(rPS(n,:))];
    
    sVp_mean = [sVp_mean; get_entry_mean(sVp(n,:))];
    rVp_mean = [rVp_mean; get_entry_mean(rVp(n,:))];
    
    sVb_mean = [sVb_mean; get_entry_mean(sVb(n,:))];
    rVb_mean = [rVb_mean; get_entry_mean(rVb(n,:))];
    
    sVisf_mean = [sVisf_mean; get_entry_mean(sVisf(n,:))];
    rVisf_mean = [rVisf_mean; get_entry_mean(rVisf(n,:))];

    sKi_MF_mean = [sKi_MF_mean; get_entry_mean(sKi_MF(n,:))];
    rKi_MF_mean = [rKi_MF_mean; get_entry_mean(rKi_MF(n,:))];
    
    sKi_Fermi_mean = [sKi_Fermi_mean; get_entry_mean(sKi_Fermi(n,:))];
    rKi_Fermi_mean = [rKi_Fermi_mean; get_entry_mean(rKi_Fermi(n,:))];
    
    sKi_TwoCompExp_mean = [sKi_TwoCompExp_mean; get_entry_mean(sKi_TwoCompExp(n,:))];
    rKi_TwoCompExp_mean = [rKi_TwoCompExp_mean; get_entry_mean(rKi_TwoCompExp(n,:))];
    
    sKi_BTEX_mean = [sKi_BTEX_mean; get_entry_mean(sKi_BTEX(n,:))];
    rKi_BTEX_mean = [rKi_BTEX_mean; get_entry_mean(rKi_BTEX(n,:))]; 
    
    sDelay_mean = [sDelay_mean; get_entry_mean(sDelay(n,:))];
    rDelay_mean = [rDelay_mean; get_entry_mean(rDelay(n,:))];
    
    sFermiDelay_mean = [sFermiDelay_mean; get_entry_mean(sFermiDelay(n,:))];
    rFermiDelay_mean = [rFermiDelay_mean; get_entry_mean(rFermiDelay(n,:))];
end

res_table = table(scanInd, patientID, scanDate, scanTime, age, gender, stressHB, restHB, hematocrit, LVMass, LVEDV, ... 
                sf, rf, sf_i, rf_i, ... 
                sE, rE, sE_i, rE_i, ... 
                sPS, rPS, sPS_i, rPS_i, ...
                sVisf, rVisf, sVisf_i, rVisf_i, ... 
                sVp, rVp, sVp_i, rVp_i, ... 
                sVb, rVb, sVb_i, rVb_i, ... 
                sKi_MF, rKi_MF, sKi_MF_i, rKi_MF_i, ... 
                sKi_Fermi, rKi_Fermi, sKi_Fermi_i, rKi_Fermi_i, ... 
                sKi_TwoCompExp, rKi_TwoCompExp, sKi_TwoCompExp_i, rKi_TwoCompExp_i, ...
                sKi_BTEX, rKi_BTEX, sKi_BTEX_i, rKi_BTEX_i, ... 
                sSD, rSD, sSD_i, rSD_i, sSplenic, rSplenic, pre_T1_blood, post_T1_blood, pre_T1_myo, post_T1_myo, ... 
                sDelay, rDelay, sFermiDelay, rFermiDelay, ... 
                ...
                sf_mean, rf_mean, sE_mean, rE_mean, sPS_mean, rPS_mean, sVp_mean, rVp_mean, sVb_mean, rVb_mean, sVisf_mean, rVisf_mean, ...
                sKi_MF_mean, rKi_MF_mean, sKi_Fermi_mean, rKi_Fermi_mean, sKi_TwoCompExp_mean, rKi_TwoCompExp_mean, sKi_BTEX_mean, rKi_BTEX_mean, ... 
                sDelay_mean, rDelay_mean, sFermiDelay_mean, rFermiDelay_mean, ...
                ...
                sSNR, rSNR, sGd, rGd, sT2S, rT2S, sAifPeakGd, rAifPeakGd, ...
                age, gender, stressHB, restHB, height, weight, diabetes, HTN, Chol, smoker, EF);

disp('=======================================================================');
disp(['Stress flow - ' num2str(mean(sf(:))) '+/-' num2str(std(sf(:)))]);
disp(['Rest flow - ' num2str(mean(rf(:))) '+/-' num2str(std(rf(:)))]);

ind = find(sf_i(:)>0);
disp(['Stress flow, ischemia - ' num2str(mean(sf_i(ind))) '+/-' num2str(std(sf_i(ind)))]);

end

function res = PerformGadgetronRecon_Matlab_ROIValues_OneCase(s1, s2, s3, a1, a2, a3, HCT)
% res = PerformGadgetronRecon_Matlab_ROIValues_OneCase(s1, s2, s3, a1, a2, a3)
% given the ROIs s1/s2/s3, get the flow and other values
% res has [flow, Ki, E, Visf, Vp, PS, SD, Ki_MF, Ki_Fermi, Ki_TwoCompExp, Ki_BTEX]

if(isfield(a1, 'flowmaps_grappa_PSIR'))
    
    [f1, f2, f3] = get_roi_values(a1.flowmaps_grappa_PSIR, a2.flowmaps_grappa_PSIR, a3.flowmaps_grappa_PSIR, s1, s2, s3);
    res.flow = [f1.m f2.m f3.m];

    EMap1 = a1.Ki_whole_grappa_PSIR(:,:,4) ./ (a1.flowmaps_grappa_PSIR(:,:,end)+eps);
    EMap2 = a2.Ki_whole_grappa_PSIR(:,:,4) ./ (a2.flowmaps_grappa_PSIR(:,:,end)+eps);
    EMap3 = a3.Ki_whole_grappa_PSIR(:,:,4) ./ (a3.flowmaps_grappa_PSIR(:,:,end)+eps);
    
    [f1, f2, f3] = get_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
    res.E = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.grappa_interVolumeMap_grappa_PSIR, a2.grappa_interVolumeMap_grappa_PSIR, a3.grappa_interVolumeMap_grappa_PSIR, s1, s2, s3);
    res.Visf = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.blood_volume_maps_grappa_PSIR, a2.blood_volume_maps_grappa_PSIR, a3.blood_volume_maps_grappa_PSIR, s1, s2, s3);
    Vp_1 = median(f1.data(:));
    Vp_2 = median(f2.data(:));
    Vp_3 = median(f3.data(:));
    res.Vp = [Vp_1 Vp_2 Vp_3] * (1-HCT);
    res.Vb = [Vp_1 Vp_2 Vp_3];

    [f1, f2, f3] = get_roi_values(a1.PS_maps_grappa_PSIR, a2.PS_maps_grappa_PSIR, a3.PS_maps_grappa_PSIR, s1, s2, s3);
    res.PS = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.SD_maps_grappa_PSIR, a2.SD_maps_grappa_PSIR, a3.SD_maps_grappa_PSIR, s1, s2, s3);
    res.SD = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,1)), s1, s2, s3);
    res.Ki_MF = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,2)), s1, s2, s3);
    res.Ki_Fermi = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,3)), s1, s2, s3);
    res.Ki_TwoCompExp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,4)), s1, s2, s3);
    res.Ki_BTEX = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( a1.btex_delay, a2.btex_delay, a3.btex_delay, s1, s2, s3);
    res.delay = [f1.m f2.m f3.m];
    
    [f1, f2, f3] = get_roi_values( a1.fermi_delay, a2.fermi_delay, a3.fermi_delay, s1, s2, s3);
    res.fermi_delay = [f1.m f2.m f3.m];
    
    %% if having second roi
    if( (~isempty(s1) & numel(s1.ROI_info_table)==2) | (~isempty(s2) & numel(s2.ROI_info_table)==2) | (~isempty(s3) & numel(s3.ROI_info_table)==2))

        [f1, f2, f3] = get_2nd_roi_values(a1.flowmaps_grappa_PSIR, a2.flowmaps_grappa_PSIR, a3.flowmaps_grappa_PSIR, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.grappa_interVolumeMap_grappa_PSIR, a2.grappa_interVolumeMap_grappa_PSIR, a3.grappa_interVolumeMap_grappa_PSIR, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.blood_volume_maps_grappa_PSIR, a2.blood_volume_maps_grappa_PSIR, a3.blood_volume_maps_grappa_PSIR, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m] * (1-HCT);
        res.Vb_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.PS_maps_grappa_PSIR, a2.PS_maps_grappa_PSIR, a3.PS_maps_grappa_PSIR, s1, s2, s3);
        res.PS_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.SD_maps_grappa_PSIR, a2.SD_maps_grappa_PSIR, a3.SD_maps_grappa_PSIR, s1, s2, s3);
        res.SD_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,1)), s1, s2, s3);
        res.Ki_MF_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,2)), s1, s2, s3);
        res.Ki_Fermi_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,3)), s1, s2, s3);
        res.Ki_TwoCompExp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,4)), s1, s2, s3);
        res.Ki_BTEX_i = [f1.m f2.m f3.m];
    end
elseif(isfield(a1, 'flowmaps_grappa_PSIR_OnlyGlobalSearch'))
    
    [f1, f2, f3] = get_roi_values(a1.flowmaps_grappa_PSIR_OnlyGlobalSearch, a2.flowmaps_grappa_PSIR_OnlyGlobalSearch, a3.flowmaps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
    res.flow = [f1.m f2.m f3.m];

    EMap1 = a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4) ./ (a1.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end)+eps);
    EMap2 = a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4) ./ (a2.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end)+eps);
    EMap3 = a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4) ./ (a3.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end)+eps);
    
    [f1, f2, f3] = get_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
    res.E = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, a2.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, a3.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
    res.Visf = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, a2.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, a3.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
    res.Vp = [f1.m f2.m f3.m] * (1-HCT);

    [f1, f2, f3] = get_roi_values(a1.PS_maps_grappa_PSIR_OnlyGlobalSearch, a2.PS_maps_grappa_PSIR_OnlyGlobalSearch, a3.PS_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
    res.PS = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.SD_maps_grappa_PSIR_OnlyGlobalSearch, a2.SD_maps_grappa_PSIR_OnlyGlobalSearch, a3.SD_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
    res.SD = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), s1, s2, s3);
    res.Ki_MF = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), s1, s2, s3);
    res.Ki_Fermi = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), s1, s2, s3);
    res.Ki_TwoCompExp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), s1, s2, s3);
    res.Ki_BTEX = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( a1.btex_delay, a2.btex_delay, a3.btex_delay, s1, s2, s3);
    res.delay = [f1.m f2.m f3.m];
    
    [f1, f2, f3] = get_roi_values( a1.fermi_delay, a2.fermi_delay, a3.fermi_delay, s1, s2, s3);
    res.fermi_delay = [f1.m f2.m f3.m];

    %% if having second roi
    if( (~isempty(s1) & numel(s1.ROI_info_table)==2) | (~isempty(s2) & numel(s2.ROI_info_table)==2) | (~isempty(s3) & numel(s3.ROI_info_table)==2))

        [f1, f2, f3] = get_2nd_roi_values(a1.flowmaps_grappa_PSIR_OnlyGlobalSearch, a2.flowmaps_grappa_PSIR_OnlyGlobalSearch, a3.flowmaps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, a2.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, a3.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, a2.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, a3.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m] * (1-HCT);

        [f1, f2, f3] = get_2nd_roi_values(a1.PS_maps_grappa_PSIR_OnlyGlobalSearch, a2.PS_maps_grappa_PSIR_OnlyGlobalSearch, a3.PS_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
        res.PS_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.SD_maps_grappa_PSIR_OnlyGlobalSearch, a2.SD_maps_grappa_PSIR_OnlyGlobalSearch, a3.SD_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
        res.SD_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), s1, s2, s3);
        res.Ki_MF_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), s1, s2, s3);
        res.Ki_Fermi_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), s1, s2, s3);
        res.Ki_TwoCompExp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), s1, s2, s3);
        res.Ki_BTEX_i = [f1.m f2.m f3.m];
    end
elseif(isfield(a1, 'flowmaps_grappa_PSIR_without_R2Star'))
    
    [f1, f2, f3] = get_roi_values(a1.flowmaps_grappa_PSIR_without_R2Star, a2.flowmaps_grappa_PSIR_without_R2Star, a3.flowmaps_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.flow = [f1.m f2.m f3.m];

    EMap1 = a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,4) ./ (a1.flowmaps_grappa_PSIR_without_R2Star(:,:,end)+eps);
    EMap2 = a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,4) ./ (a2.flowmaps_grappa_PSIR_without_R2Star(:,:,end)+eps);
    EMap3 = a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,4) ./ (a3.flowmaps_grappa_PSIR_without_R2Star(:,:,end)+eps);
    
    [f1, f2, f3] = get_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
    res.E = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a2.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a3.grappa_interVolumeMap_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.Visf = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.blood_volume_maps_grappa_PSIR_without_R2Star, a2.blood_volume_maps_grappa_PSIR_without_R2Star, a3.blood_volume_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.Vp = [f1.m f2.m f3.m] * (1-HCT);
    res.Vb = [f1.m f2.m f3.m];
    
    [f1, f2, f3] = get_roi_values(a1.PS_maps_grappa_PSIR_without_R2Star, a2.PS_maps_grappa_PSIR_without_R2Star, a3.PS_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.PS = [f1.m f2.m f3.m];

    if(isfield(a1, 'SD_maps_grappa_PSIR_without_R2Star'))
        [f1, f2, f3] = get_roi_values(a1.SD_maps_grappa_PSIR_without_R2Star, a2.SD_maps_grappa_PSIR_without_R2Star, a3.SD_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.SD = [f1.m f2.m f3.m];
    else
        res.SD = [-1 -1 -1];
    end
    
    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), s1, s2, s3);
    res.Ki_MF = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), s1, s2, s3);
    res.Ki_Fermi = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), s1, s2, s3);
    res.Ki_TwoCompExp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), s1, s2, s3);
    res.Ki_BTEX = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( a1.btex_delay, a2.btex_delay, a3.btex_delay, s1, s2, s3);
    res.delay = [f1.m f2.m f3.m];
    
    [f1, f2, f3] = get_roi_values( a1.fermi_delay, a2.fermi_delay, a3.fermi_delay, s1, s2, s3);
    res.fermi_delay = [f1.m f2.m f3.m];

    %% if having second roi
    if( (~isempty(s1) & numel(s1.ROI_info_table)==2) | (~isempty(s2) & numel(s2.ROI_info_table)==2) | (~isempty(s3) & numel(s3.ROI_info_table)==2))

        [f1, f2, f3] = get_2nd_roi_values(a1.flowmaps_grappa_PSIR_without_R2Star, a2.flowmaps_grappa_PSIR_without_R2Star, a3.flowmaps_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a2.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a3.grappa_interVolumeMap_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.blood_volume_maps_grappa_PSIR_without_R2Star, a2.blood_volume_maps_grappa_PSIR_without_R2Star, a3.blood_volume_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m] * (1-HCT);
        res.Vb_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.PS_maps_grappa_PSIR_without_R2Star, a2.PS_maps_grappa_PSIR_without_R2Star, a3.PS_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.PS_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.SD_maps_grappa_PSIR_without_R2Star, a2.SD_maps_grappa_PSIR_without_R2Star, a3.SD_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.SD_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), s1, s2, s3);
        res.Ki_MF_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), s1, s2, s3);
        res.Ki_Fermi_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), s1, s2, s3);
        res.Ki_TwoCompExp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), s1, s2, s3);
        res.Ki_BTEX_i = [f1.m f2.m f3.m];
    end
end

end

function [f1, f2, f3] = get_roi_values(a1, a2, a3, s1, s2, s3)
    if(isempty(s1) | isempty(a1))
        f1.m = -1;
        f1.data = [-1];
    else
        f1 = roi_statistics(flipdim(a1(:,:,end), 2), s1.ROI_info_table(1,1));
    end
    
    if(isempty(s2) | isempty(a2))
        f2.m = -1;
        f2.data = [-1];
    else
        f2 = roi_statistics(flipdim(a2(:,:,end), 2), s2.ROI_info_table(1,1));
    end
    
    if(isempty(s3) | isempty(a3))
        f3.m = -1;
        f3.data = [-1];
    else
        f3 = roi_statistics(flipdim(a3(:,:,end), 2), s3.ROI_info_table(1,1));
    end
end

function [f1, f2, f3] = get_2nd_roi_values(a1, a2, a3, s1, s2, s3)

    if(~isempty(s1) & numel(s1.ROI_info_table)==2)
        f1 = roi_statistics(flipdim(a1(:,:,end), 2), s1.ROI_info_table(2,1));
    else
        f1.m = -1;
    end
    
    if(~isempty(s2) & numel(s2.ROI_info_table)==2)
        f2 = roi_statistics(flipdim(a2(:,:,end), 2), s2.ROI_info_table(2,1));
    else
        f2.m = -1;
    end
    
    if(~isempty(s3) & numel(s3.ROI_info_table)==2)
        f3 = roi_statistics(flipdim(a3(:,:,end), 2), s3.ROI_info_table(2,1));
    else
        f3.m = -1;
    end
end

function v = get_valid_value(v_in)

    if(isnan(v_in))
        v = 0;
    else
        if(isstr(v_in))
            v = 1;
        else
            v = v_in;
        end
    end
end

function v = get_entry_mean(f)
    ind = find(f>0);
    if(isempty(ind))
        v=-1;
    else
        v = mean(f(ind));
    end
end
