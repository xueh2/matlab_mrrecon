
function [ROITable, sf, rf, sf_i, rf_i, res_table] = PerformGadgetronRecon_SavedIsmrmrd_Matlab_PlotPerf(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, fixed_HCT, prefix, processing_always)
% [ROITable, sf, rf, sf_i, rf_i, res_table] = PerformGadgetronRecon_SavedIsmrmrd_Matlab_ROIValues(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, prefix, processing_always)
% if isempty(fixed_HCT)==1, measured HCT is used; if(fixed_HCT>=0 & fixed_HCT<=1), this value is used

if(nargin<9)
    processing_always = 1;
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
stressHB = [];
restHB = [];
hematocrit = [];

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
    
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]); 
    disp(['==================================================================']);  
       
    [configName, scannerID, pID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(stressCase);

    figDir = fullfile(resDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_Figure']) 
    
    roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_ROI'])
    
    scanInd = [scanInd; n];
    patientID = [patientID; {pID}];
    scanDate = [scanDate; study_dates];
    scanTime = [scanTime; study_time];
    age = [age; PerfTable{n+1, 14}];
    gender = [gender; PerfTable{n+1, 15}];
    stressHB = [stressHB; PerfTable{n+1, 9}];
    restHB = [restHB; PerfTable{n+1, 13}];
    HCT = PerfTable{n+1, hct_column};
    if(isnan(HCT))
        HCT = 0;
    else
        if(HCT>1)
            HCT = HCT/100;
        end
    end
    hematocrit = [hematocrit; HCT];
    
    if(~isempty(fixed_HCT))
        if(fixed_HCT>=0 & fixed_HCT<=1)
            HCT = fixed_HCT;
        end
    end
    
    disp([num2str(n) ' out of ' num2str(num) ' - HCT : ' num2str(HCT)]); 
   
    hct_str = num2str(HCT);
    ind = find(hct_str=='.');
    if(~isempty(ind))
        hct_str(ind(:)) = 'p';
    end
    suffix = [study_dates '_hct' hct_str];
    
    if(isFileExist(fullfile(roiDir, 'myo_stress1.mat')))
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
        
    if(isFileExist(fullfile(roiDir, s1_roi)) & isFileExist(fullfile(roiDir, r1_roi)))
        disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]);       
        
        two_ROI = 0;

        PerfResult_file = fullfile(roiDir, [prefix '_PerfResult_' suffix '.mat']);
        has_result_file = 0;
        if(~processing_always & isFileExist(PerfResult_file))
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
            
            if(two_ROI)
                sE_i = [sE_i; res_stress.E_i];
                rE_i = [rE_i; res_rest.E_i];

                sPS_i = [sPS_i; res_stress.PS_i];
                rPS_i = [rPS_i; res_rest.PS_i];

                sVisf_i = [sVisf_i; res_stress.Visf_i];
                rVisf_i = [rVisf_i; res_rest.Visf_i];

                sVp_i = [sVp_i; res_stress.Vp_i];
                rVp_i = [rVp_i; res_rest.Vp_i];

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
        else        
            cd(roiDir)
            v = PerfTable(n+1, :);

            nV = numel(v);

            s1 = load(fullfile(roiDir, s1_roi));
            s2 = load(fullfile(roiDir, s2_roi));
            s3 = load(fullfile(roiDir, s3_roi));

            r1 = load(fullfile(roiDir, r1_roi));
            r2 = load(fullfile(roiDir, r2_roi));
            r3 = load(fullfile(roiDir, r3_roi));

            if(~isFileExist(fullfile(roiDir, 'rest.mat')))
                copyfile(fullfile(figDir, 'rest.mat'), roiDir);
            end
            
            if(~isFileExist(fullfile(roiDir, 'stress.mat')))
                copyfile(fullfile(figDir, 'stress.mat'), roiDir);
            end
            
            rest = load(fullfile(roiDir, 'rest.mat'));
            stress = load(fullfile(roiDir, 'stress.mat'));
            
            stressMat1 = load(fullfile(resDir, study_dates, stressCase, [prefix '_' suffix '_0.mat']));
            stressMat2 = load(fullfile(resDir, study_dates, stressCase, [prefix '_' suffix '_1.mat']));
            stressMat3 = load(fullfile(resDir, study_dates, stressCase, [prefix '_' suffix '_2.mat']));
            
            restMat1 = load(fullfile(resDir, study_dates, restCase, [prefix '_' suffix '_0.mat']));
            restMat2 = load(fullfile(resDir, study_dates, restCase, [prefix '_' suffix '_1.mat']));
            restMat3 = load(fullfile(resDir, study_dates, restCase, [prefix '_' suffix '_2.mat']));
                  
            % figure; imagescn(stress.flow_stress, [0 6]); PerfColorMap; figure; imagescn(stressMat1.flowmaps_grappa_PSIR(:,:,end), [0 6]); PerfColorMap; figure; imagescn(stressMat1.Ki_whole_grappa_PSIR, [0 6]); PerfColorMap;
            % figure; imagescn(stress.flow_stress, [0 6]); PerfColorMap; figure; imagescn(stressMat3.flowmaps_grappa_PSIR(:,:,end), [0 6]); PerfColorMap; figure; imagescn(stressMat3.Ki_whole_grappa_PSIR, [0 6]); PerfColorMap;

            res_stress = PerformGadgetronRecon_Matlab_ROIValues_OneCase(s1, s2, s3, stressMat1, stressMat2, stressMat3, HCT);
            res_rest = PerformGadgetronRecon_Matlab_ROIValues_OneCase(r1, r2, r3, restMat1, restMat2, restMat3, HCT);
            
            v{nV+1} = res_stress.flow;
            disp(['stress flow : ' num2str(v{nV+1})]);

            v{nV+2} = res_rest.flow;
            disp(['rest flow : ' num2str(v{nV+2})]);

            if(~isempty(find(res_stress.flow>5.5)))
                figure; imagescn(stressMat3.flowmaps_grappa_PSIR(:,:,end), [0 8]); PerfColorMap;
                figure; imagescn(stressMat2.flowmaps_grappa_PSIR(:,:,end), [0 8]); PerfColorMap;
                figure; imagescn(stressMat3.flowmaps_grappa_PSIR(:,:,end), [0 8]); PerfColorMap;
                pause;
                closeall
            end
            
            if(numel(s1.ROI_info_table)==2)
                
                two_ROI = 1;

                v{nV+3} = res_stress.flow_i;
                disp(['stress flow, ischemia : ' num2str(v{nV+3})]);
                
                v{nV+4} = res_rest.flow_i;
            else
                v{nV+3} = [-1 -1 -1];
                v{nV+4} = [-1 -1 -1];
            end

            v{nV+5} = max(stress.aif_stress_baseline_corrected);
            v{nV+6} = max(rest.aif_rest_baseline_corrected);
            v{nV+7} = HCT;
            
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

            if(two_ROI)
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
            save(PerfResult_file, 'PerfResult', 'res_stress', 'res_rest', 'two_ROI');           
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

res_table = table(scanInd, patientID, scanDate, scanTime, age, gender, stressHB, restHB, hematocrit, ... 
                sf, rf, sf_i, rf_i, ... 
                sE, rE, sE_i, rE_i, ... 
                sPS, rPS, sPS_i, rPS_i, ...
                sVisf, rVisf, sVisf_i, rVisf_i, ... 
                sVp, rVp, sVp_i, rVp_i, ... 
                sKi_MF, rKi_MF, sKi_MF_i, rKi_MF_i, ... 
                sKi_Fermi, rKi_Fermi, sKi_Fermi_i, rKi_Fermi_i, ... 
                sKi_TwoCompExp, rKi_TwoCompExp, sKi_TwoCompExp_i, rKi_TwoCompExp_i, ...
                sKi_BTEX, rKi_BTEX, sKi_BTEX_i, rKi_BTEX_i, ... 
                sSD, rSD, sSD_i, rSD_i);

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
    res.Vp = [f1.m f2.m f3.m] * (1-HCT);

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

    %% if having second roi
    if(numel(s1.ROI_info_table)==2)

        [f1, f2, f3] = get_2nd_roi_values(a1.flowmaps_grappa_PSIR, a2.flowmaps_grappa_PSIR, a3.flowmaps_grappa_PSIR, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.E_stress, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.grappa_interVolumeMap_grappa_PSIR, a2.grappa_interVolumeMap_grappa_PSIR, a3.grappa_interVolumeMap_grappa_PSIR, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.blood_volume_maps_grappa_PSIR, a2.blood_volume_maps_grappa_PSIR, a3.blood_volume_maps_grappa_PSIR, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m] * (1-HCT);

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
    
    [f1, f2, f3] = get_roi_values(stress.E_stress, s1, s2, s3);
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

    %% if having second roi
    if(numel(s1.ROI_info_table)==2)

        [f1, f2, f3] = get_2nd_roi_values(a1.flowmaps_grappa_PSIR_OnlyGlobalSearch, a2.flowmaps_grappa_PSIR_OnlyGlobalSearch, a3.flowmaps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.E_stress, s1, s2, s3);
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
    
    [f1, f2, f3] = get_roi_values(stress.E_stress, s1, s2, s3);
    res.E = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a2.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a3.grappa_interVolumeMap_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.Visf = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.blood_volume_maps_grappa_PSIR_without_R2Star, a2.blood_volume_maps_grappa_PSIR_without_R2Star, a3.blood_volume_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.Vp = [f1.m f2.m f3.m] * (1-HCT);

    [f1, f2, f3] = get_roi_values(a1.PS_maps_grappa_PSIR_without_R2Star, a2.PS_maps_grappa_PSIR_without_R2Star, a3.PS_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.PS = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.SD_maps_grappa_PSIR_without_R2Star, a2.SD_maps_grappa_PSIR_without_R2Star, a3.SD_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.SD = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), s1, s2, s3);
    res.Ki_MF = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), s1, s2, s3);
    res.Ki_Fermi = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), s1, s2, s3);
    res.Ki_TwoCompExp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), s1, s2, s3);
    res.Ki_BTEX = [f1.m f2.m f3.m];

    %% if having second roi
    if(numel(s1.ROI_info_table)==2)

        [f1, f2, f3] = get_2nd_roi_values(a1.flowmaps_grappa_PSIR_without_R2Star, a2.flowmaps_grappa_PSIR_without_R2Star, a3.flowmaps_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.E_stress, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a2.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a3.grappa_interVolumeMap_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.blood_volume_maps_grappa_PSIR_without_R2Star, a2.blood_volume_maps_grappa_PSIR_without_R2Star, a3.blood_volume_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m] * (1-HCT);

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
    f1 = roi_statistics(flipdim(a1(:,:,end), 2), s1.ROI_info_table(1,1));
    f2 = roi_statistics(flipdim(a2(:,:,end), 2), s2.ROI_info_table(1,1));
    f3 = roi_statistics(flipdim(a3(:,:,end), 2), s3.ROI_info_table(1,1));
end

function [f1, f2, f3] = get_2nd_roi_values(a1, a2, a3, s1, s2, s3)

    if(numel(s1.ROI_info_table)==2)
        f1 = roi_statistics(flipdim(a1(:,:,end), 2), s1.ROI_info_table(1,1));
    else
        f1.m = -1;
    end
    
    if(numel(s2.ROI_info_table)==2)
        f2 = roi_statistics(flipdim(a2(:,:,end), 2), s2.ROI_info_table(1,1));
    else
        f2.m = -1;
    end
    
    if(numel(s3.ROI_info_table)==2)
        f3 = roi_statistics(flipdim(a3(:,:,end), 2), s3.ROI_info_table(1,1));
    else
        f3.m = -1;
    end
end

