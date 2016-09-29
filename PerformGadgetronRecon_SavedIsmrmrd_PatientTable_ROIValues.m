
function [ROITable, sf, rf, sf_i, rf_i] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, excluded, processing_always)
% [ROITable, sf, rf, sf_i, rf_i] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, excluded, processing_always)
% [ROITable, sf, rf, sf_i, rf_i] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area', 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area_ROI')
% [ROITable, sf, rf, sf_i, rf_i] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area_Native_T1_0', 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area_ROI')
% [ROITable, sf, rf, sf_i, rf_i] = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ROIValues(PerfTable, 'I:\ReconResults\BARTS', 'D:\data\ut\NewData\PaperResults\Barts_ROI')

if(nargin<7)
    excluded = 0;
end

if(nargin<8)
    processing_always = 1;
end

ROITable = PerfTable(1, :);
ROITable = [ROITable 'stress flow' 'rest flow' 'stress ischemia flow' 'rest ischemia flow' 'stress aif peak Gd with baseline correction' 'rest aif peak Gd with baseline correction'];

nV = numel(PerfTable(1, :));

sf = [];
rf = [];
sf_i = [];
rf_i = [];

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
       
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(stressCase);

    figDir = fullfile(resDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' patientID '_' studyID '_' study_dates '_Figure']) 
    
    if(excluded)
        roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' patientID '_' studyID '_' study_dates '_ROI_excluded'])
    else
        roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' patientID '_' studyID '_' study_dates '_ROI'])
    end
    
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

        PerfResult_file = fullfile(roiDir, 'PerfResult.mat');
        if(~processing_always & isFileExist(PerfResult_file))
            tt = load(PerfResult_file);
            v = tt.PerfResult;
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

            f1 = roi_statistics(stress.flow_stress(:,:,1), s1.ROI_info_table(1,1));
            f2 = roi_statistics(stress.flow_stress(:,:,2), s2.ROI_info_table(1,1));
            f3 = roi_statistics(stress.flow_stress(:,:,3), s3.ROI_info_table(1,1));

            v{nV+1} = [f1.m f2.m f3.m];
            disp(['stress flow : ' num2str(v{nV+1})]);

            f1 = roi_statistics(rest.flow_rest(:,:,1), r1.ROI_info_table(1,1));
            f2 = roi_statistics(rest.flow_rest(:,:,2), r2.ROI_info_table(1,1));
            f3 = roi_statistics(rest.flow_rest(:,:,3), r3.ROI_info_table(1,1));

            v{nV+2} = [f1.m f2.m f3.m];
            disp(['rest flow : ' num2str(v{nV+2})]);

            if(numel(s1.ROI_info_table)==2)
                
                two_ROI = 1;
                
                if(numel(s1.ROI_info_table)==2)
                    f1 = roi_statistics(stress.flow_stress(:,:,1), s1.ROI_info_table(2));
                else
                    f1.m = -1;
                end

                if(numel(s2.ROI_info_table)==2)
                    f2 = roi_statistics(stress.flow_stress(:,:,2), s2.ROI_info_table(2));
                else
                    f2.m = -1;
                end

                if(numel(s3.ROI_info_table)==2)
                    f3 = roi_statistics(stress.flow_stress(:,:,3), s3.ROI_info_table(2));
                else
                    f3.m = -1;
                end

                v{nV+3} = [f1.m f2.m f3.m];
                disp(['stress flow, ischemia : ' num2str(v{nV+3})]);
            else
                v{nV+3} = [];
                v{nV+4} = [];
            end

            v{nV+5} = max(stress.aif_stress_baseline_corrected);
            v{nV+6} = max(rest.aif_rest_baseline_corrected);
        end
    
        ROITable = [ROITable; v];
                
        rf = [rf ROITable{end, nV+2}];       
        if(two_ROI)
            sf_i = [sf_i ROITable{end, nV+3}];
        end
        sf = [sf ROITable{end, nV+1}];
        
        PerfResult = v;
        save(fullfile(roiDir, 'PerfResult.mat'), 'PerfResult');
    else
        disp(['Cannot find the ROI ... ']);
        winopen(figDir);
        mkdir(roiDir);
        winopen(roiDir);
        pause;
    end
end

disp('=======================================================================');
disp(['Stress flow - ' num2str(mean(sf)) '+/-' num2str(std(sf))]);
disp(['Rest flow - ' num2str(mean(rf)) '+/-' num2str(std(rf))]);

ind = find(sf_i>0);
disp(['Stress flow, ischemia - ' num2str(mean(sf_i(ind))) '+/-' num2str(std(sf_i(ind)))]);

