
function PerformGadgetronRecon_SavedIsmrmrd_EndoEpi_AHAModel(PerfTable, resDir, contourDir, stress_column, rest_column)
% PerformGadgetronRecon_SavedIsmrmrd_EndoEpi_AHAModel(PerfTable, resDir, contourDir, stress_column, rest_column)

if(exist(contourDir)~=7)
    mkdir(contourDir);
end

num = size(PerfTable, 1)-1;
for n=1:num
    disp(['=============================================================================================================================================']);  
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]); 
    PerfTable{n+1, :}
    
    stressCase = PerfTable{n+1, stress_column};
    restCase = PerfTable{n+1, rest_column};
    
    process = 1;
       
    if(~process)
        continue;
    end
            
    [configName, scannerID, pID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(stressCase);
    [configName, scannerID, pID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(restCase);
    figDir = fullfile(resDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_' study_time_stress '_' study_time_rest '_Figure']);
    
    roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_ROI'])
    
    if(isFileExist(roiDir) & isFileExist(fullfile(roiDir, 'r1.mat')) ... 
            & isFileExist(fullfile(roiDir, 'r2.mat')) ...
            & isFileExist(fullfile(roiDir, 'r3.mat')) ...
            & isFileExist(fullfile(roiDir, 's1.mat')) ...
            & isFileExist(fullfile(roiDir, 's2.mat')) ... 
            & isFileExist(fullfile(roiDir, 's3.mat')) ...
            & isFileExist(fullfile(roiDir, 'r1_p.mat')) ...
            & isFileExist(fullfile(roiDir, 'r2_p.mat')) ...
            & isFileExist(fullfile(roiDir, 'r3_p.mat')) ...
            & isFileExist(fullfile(roiDir, 's1_p.mat')) ...
            & isFileExist(fullfile(roiDir, 's2_p.mat')) ...
            & isFileExist(fullfile(roiDir, 's3_p.mat')) ...
        )
        dir(roiDir);
        onlyReview = 1;
%         continue; % temporary
    else
        ls(roiDir)
        continue;
    end
    cd(roiDir)
    
    r1 = load('r1');
    r2 = load('r2');
    r3 = load('r3');
    r1_p = load('r1_p');
    r2_p = load('r2_p');
    r3_p = load('r3_p');
    
    s1 = load('s1');
    s2 = load('s2');
    s3 = load('s3');
    s1_p = load('s1_p');
    s2_p = load('s2_p');
    s3_p = load('s3_p');
    
    % ---------------------------------------
    
    endo_r1 = [r1.ROI_info_table(1).ROI_x_coordinates(:) r1.ROI_info_table(1).ROI_y_coordinates(:)];
    epi_r1 = [r1.ROI_info_table(2).ROI_x_coordinates(:) r1.ROI_info_table(2).ROI_y_coordinates(:)];
    lv_center_r1 = [r1_p.Profile_info_table(1).Profile_x_coordinates r1_p.Profile_info_table(1).Profile_y_coordinates];
    rv_insertion_r1 = [r1_p.Profile_info_table(2).Profile_x_coordinates r1_p.Profile_info_table(2).Profile_y_coordinates];
    
    endo_s1 = [s1.ROI_info_table(1).ROI_x_coordinates(:) s1.ROI_info_table(1).ROI_y_coordinates(:)];
    epi_s1 = [s1.ROI_info_table(2).ROI_x_coordinates(:) s1.ROI_info_table(2).ROI_y_coordinates(:)];
    lv_center_s1 = [s1_p.Profile_info_table(1).Profile_x_coordinates s1_p.Profile_info_table(1).Profile_y_coordinates];
    rv_insertion_s1 = [s1_p.Profile_info_table(2).Profile_x_coordinates s1_p.Profile_info_table(2).Profile_y_coordinates];
    
    % ---------------------------------------
    
    endo_r2 = [r2.ROI_info_table(1).ROI_x_coordinates(:) r2.ROI_info_table(1).ROI_y_coordinates(:)];
    epi_r2 = [r2.ROI_info_table(2).ROI_x_coordinates(:) r2.ROI_info_table(2).ROI_y_coordinates(:)];
    lv_center_r2 = [r2_p.Profile_info_table(1).Profile_x_coordinates r2_p.Profile_info_table(1).Profile_y_coordinates];
    rv_insertion_r2 = [r2_p.Profile_info_table(2).Profile_x_coordinates r2_p.Profile_info_table(2).Profile_y_coordinates];
    
    endo_s2 = [s2.ROI_info_table(1).ROI_x_coordinates(:) s2.ROI_info_table(1).ROI_y_coordinates(:)];
    epi_s2 = [s2.ROI_info_table(2).ROI_x_coordinates(:) s2.ROI_info_table(2).ROI_y_coordinates(:)];
    lv_center_s2 = [s2_p.Profile_info_table(1).Profile_x_coordinates s2_p.Profile_info_table(1).Profile_y_coordinates];
    rv_insertion_s2 = [s2_p.Profile_info_table(2).Profile_x_coordinates s2_p.Profile_info_table(2).Profile_y_coordinates];
    
    % ---------------------------------------
    
    endo_r3 = [r3.ROI_info_table(1).ROI_x_coordinates(:) r3.ROI_info_table(1).ROI_y_coordinates(:)];
    epi_r3 = [r3.ROI_info_table(2).ROI_x_coordinates(:) r3.ROI_info_table(2).ROI_y_coordinates(:)];
    lv_center_r3 = [r3_p.Profile_info_table(1).Profile_x_coordinates r3_p.Profile_info_table(1).Profile_y_coordinates];
    rv_insertion_r3 = [r3_p.Profile_info_table(2).Profile_x_coordinates r3_p.Profile_info_table(2).Profile_y_coordinates];
    
    endo_s3 = [s3.ROI_info_table(1).ROI_x_coordinates(:) s3.ROI_info_table(1).ROI_y_coordinates(:)];
    epi_s3 = [s3.ROI_info_table(2).ROI_x_coordinates(:) s3.ROI_info_table(2).ROI_y_coordinates(:)];
    lv_center_s3 = [s3_p.Profile_info_table(1).Profile_x_coordinates s3_p.Profile_info_table(1).Profile_y_coordinates];
    rv_insertion_s3 = [s3_p.Profile_info_table(2).Profile_x_coordinates s3_p.Profile_info_table(2).Profile_y_coordinates];
    
    % ---------------------------------------
    
    figure; 
    subplot(2, 3, 1);
    hold on; plot(endo_s1(:,1), endo_s1(:,2), 'r');plot(epi_s1(:,1), epi_s1(:,2), 'b'); plot(lv_center_s1(1), lv_center_s1(2), 'r.'); plot(rv_insertion_s1(1), rv_insertion_s1(2), 'bs'); title('rest 1')
    subplot(2, 3, 2);
    hold on; plot(endo_r1(:,1), endo_r1(:,2), 'r');plot(epi_r1(:,1), epi_r1(:,2), 'b'); plot(lv_center_r1(1), lv_center_r1(2), 'r.'); plot(rv_insertion_r1(1), rv_insertion_r1(2), 'bs'); title('stress 1')
    subplot(2, 3, 3);
    hold on; plot(endo_s2(:,1), endo_s2(:,2), 'r');plot(epi_s2(:,1), epi_s2(:,2), 'b'); plot(lv_center_s2(1), lv_center_s2(2), 'r.'); plot(rv_insertion_s2(1), rv_insertion_s2(2), 'bs'); title('rest 2')
    subplot(2, 3, 4);
    hold on; plot(endo_r2(:,1), endo_r2(:,2), 'r');plot(epi_r2(:,1), epi_r2(:,2), 'b'); plot(lv_center_r2(1), lv_center_r2(2), 'r.'); plot(rv_insertion_r2(1), rv_insertion_r2(2), 'bs'); title('stress 2')
    subplot(2, 3, 5);
    hold on; plot(endo_s3(:,1), endo_s3(:,2), 'r');plot(epi_s3(:,1), epi_s3(:,2), 'b'); plot(lv_center_s3(1), lv_center_s3(2), 'r.'); plot(rv_insertion_s3(1), rv_insertion_s3(2), 'bs'); title('rest 3')
    subplot(2, 3, 6);
    hold on; plot(endo_r3(:,1), endo_r3(:,2), 'r');plot(epi_r3(:,1), epi_r3(:,2), 'b'); plot(lv_center_r3(1), lv_center_r3(2), 'r.'); plot(rv_insertion_r3(1), rv_insertion_r3(2), 'bs'); title('stress 3')
    
    AHA_r1 = SplitEndoEpiContourForAHAModel(endo_r1, epi_r1, lv_center_r1, rv_insertion_r1, 6);
    AHA_r2 = SplitEndoEpiContourForAHAModel(endo_r2, epi_r2, lv_center_r2, rv_insertion_r2, 6);
    AHA_r3 = SplitEndoEpiContourForAHAModel(endo_r3, epi_r3, lv_center_r3, rv_insertion_r3, 4);

    AHA_s1 = SplitEndoEpiContourForAHAModel(endo_s1, epi_s1, lv_center_s1, rv_insertion_s1, 6);
    AHA_s2 = SplitEndoEpiContourForAHAModel(endo_s2, epi_s2, lv_center_s2, rv_insertion_s2, 6);
    AHA_s3 = SplitEndoEpiContourForAHAModel(endo_s3, epi_s3, lv_center_s3, rv_insertion_s3, 4);
    
    save('AHA_model', 'AHA_r1', 'AHA_r2', 'AHA_r3', 'AHA_s1', 'AHA_s2', 'AHA_s3');
    
    figure; 
    h = subplot(2, 3, 1);
    h = PlotContourForAHAModel(AHA_r1, h); title('rest 1')
    h = subplot(2, 3, 4);
    h = PlotContourForAHAModel(AHA_s1, h); title('stress 1')
    h = subplot(2, 3, 2);
    h = PlotContourForAHAModel(AHA_r2, h); title('rest 2')
    h = subplot(2, 3, 5);
    h = PlotContourForAHAModel(AHA_s2, h); title('stress 2')
    h = subplot(2, 3, 3);
    h = PlotContourForAHAModel(AHA_r3, h); title('rest 3')
    h = subplot(2, 3, 6);
    h = PlotContourForAHAModel(AHA_s3, h); title('stress 3') 
   
    pause;
    closeall;
end
