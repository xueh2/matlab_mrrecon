
function PerformGadgetronRecon_SavedIsmrmrd_PerfCases_ContourDrawing(perf_cases, resDir, contourDir)
% PerformGadgetronRecon_SavedIsmrmrd_PerfCases_ContourDrawing(perf_cases, resDir, contourDir)

if(exist(contourDir)~=7)
    mkdir(contourDir);
end

num = size(perf_cases, 1);
start = 1;
stress_column = 2;
rest_column = 3;

for n=start:num
    disp(['=============================================================================================================================================']);  
    
    stressCase = perf_cases{n, stress_column}
    restCase = perf_cases{n, rest_column};

    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' stressCase ' - ' restCase]); 
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates_stress, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(stressCase);
    [configName, scannerID, patientID, studyID, measurementID, study_dates_rest, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(restCase);
    
    roiDir = fullfile(contourDir, study_dates_stress, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' patientID '_' studyID '_' study_dates_stress '_ROI'])
    
    if(isFileExist(roiDir) & isFileExist(fullfile(roiDir, 'r1.mat')) ... 
            & isFileExist(fullfile(roiDir, 'perf_mask_stress_0.mat')) ... 
            & isFileExist(fullfile(roiDir, 'perf_mask_stress_1.mat')) ... 
            & isFileExist(fullfile(roiDir, 'perf_mask_stress_2.mat')) ... 
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
        disp('All files are found ... ');
        
        resCaseDir = [scannerID '_' patientID '_' studyID '_' study_dates_stress];    

        stressDir = fullfile(resDir, study_dates_stress, stressCase);
        
        ind = find(stressDir=='\');
        ind2 = strfind(stressDir, scannerID);
        data_name = stressDir(ind(end)+1:ind2(1)-2);
        figDir = fullfile(resDir, study_dates_stress, [data_name '_' resCaseDir '_' study_time_stress '_' study_time_rest '_Figure']);
        disp(figDir)

        rest = load(fullfile(figDir, 'rest.mat'));
        stress = load(fullfile(figDir, 'stress.mat'));
        
        redo_flag = 0;
        
        scrsz = get(0, 'ScreenSize');
        bottom = 64;
        
        for kk=1:3
            figure('Position', [(kk-1)*512+128 bottom 512 512]); imagescn(stress.flow_stress(:,:,kk), [0 6], [], [], [], fullfile(roiDir, ['s' num2str(kk) '.mat'])); PerfColorMap;
        end
        
        try
            for kk=1:3
                figure('Position', [(kk-1)*512+128 bottom+512 512 512]); 
                imagescn(stress.flow_stress(:,:,kk), [0 6], [], [], []); PerfColorMap;
                pt = load(fullfile(roiDir, ['s' num2str(kk) '_p.mat']));
                hold on
                plot(pt.Profile_info_table(1).Profile_x_coordinates, pt.Profile_info_table(1).Profile_y_coordinates, 'r+', 'MarkerSize', 12, 'LineWidth', 4);
                plot(pt.Profile_info_table(2).Profile_x_coordinates, pt.Profile_info_table(2).Profile_y_coordinates, 'b+', 'MarkerSize', 12, 'LineWidth', 4);
                plot(pt.Profile_info_table(3).Profile_x_coordinates, pt.Profile_info_table(3).Profile_y_coordinates, 'g+', 'MarkerSize', 12, 'LineWidth', 4);
                hold off
            end
        catch
            redo_flag = 1;
        end
        
        for kk=1:3
            figure('Position', [(kk-1)*512+128 bottom 512 512]); imagescn(rest.flow_rest(:,:,kk), [0 6], [], [], [], fullfile(roiDir, ['r' num2str(kk) '.mat'])); PerfColorMap;
        end
        
        try
            for kk=1:3
                figure('Position', [(kk-1)*512+128 bottom+512 512 512]); 
                imagescn(rest.flow_rest(:,:,kk), [0 6], [], [], []); PerfColorMap;
                pt = load(fullfile(roiDir, ['r' num2str(kk) '_p.mat']));
                hold on
                plot(pt.Profile_info_table(1).Profile_x_coordinates, pt.Profile_info_table(1).Profile_y_coordinates, 'r+', 'MarkerSize', 12, 'LineWidth', 4);
                plot(pt.Profile_info_table(2).Profile_x_coordinates, pt.Profile_info_table(2).Profile_y_coordinates, 'b+', 'MarkerSize', 12, 'LineWidth', 4);
                plot(pt.Profile_info_table(3).Profile_x_coordinates, pt.Profile_info_table(3).Profile_y_coordinates, 'g+', 'MarkerSize', 12, 'LineWidth', 4);
                hold off
            end
        catch
            redo_flag = 1;
        end
        
        if(~redo_flag)
            reply = input('Do you want to redo contouring? Y/N [N]:','s');
            if isempty(reply)
               reply = 'N';
            end

            if(reply=='N')
                closeall
                continue; % temporary
            else
                onlyReview = 1;
            end        
        else
            onlyReview = 1;
        end
        
        closeall
    else    
        mkdir(roiDir)
        onlyReview = 1;
    end
    cd(roiDir)
    
    try
        [h_flow_stress, h_flow_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir,  stressCase, restCase, [0 6], onlyReview, 0);            
    catch
        continue;
    end
    
    figure;
    for s=1:numel(h_flow_stress)
        mask_file = fullfile(roiDir, ['perf_mask_stress_' num2str(s-1) '.mat']);
        if(~isFileExist(mask_file))
            [I, xlim, ylim, clim, position]=getimagedata(h_flow_stress(s));        
            BW = roipoly(I);                
            save(mask_file, 'BW');
            
            BW_stress(:,:,s) = BW;
        else
            v = load(mask_file);
            BW_stress(:,:,s) = v.BW;
        end
    end
            
   figure; imagescn(BW_stress);
   reply = input('Do you want to redo stress mask? Y/N [N]:','s');
   if isempty(reply)
      reply = 'N';
   end
   
   if(reply=='Y')
       figure
       for s=1:numel(h_flow_stress)
          mask_file = fullfile(roiDir, ['perf_mask_stress_' num2str(s-1) '.mat']);
            [I, xlim, ylim, clim, position]=getimagedata(h_flow_stress(s));        
            BW = roipoly(I);                
            save(mask_file, 'BW');

            BW_stress(:,:,s) = BW;
       end
   end
   
   figure;
    for s=1:numel(h_flow_rest)
        mask_file = fullfile(roiDir, ['perf_mask_rest_' num2str(s-1) '.mat']);
        if(~isFileExist(mask_file))
            [I, xlim, ylim, clim, position]=getimagedata(h_flow_rest(s));        
            BW = roipoly(I);                
            save(mask_file, 'BW');
            BW_rest(:,:,s) = BW;         
        else
            v = load(mask_file);
            BW_rest(:,:,s) = v.BW;            
        end
    end
    
   figure; imagescn(BW_rest);
   reply = input('Do you want to redo rest mask? Y/N [N]:','s');
   if isempty(reply)
      reply = 'N';
   end
   
   if(reply=='Y')
    figure;
    for s=1:numel(h_flow_rest)
        mask_file = fullfile(roiDir, ['perf_mask_rest_' num2str(s-1) '.mat']);
        [I, xlim, ylim, clim, position]=getimagedata(h_flow_rest(s));        
        BW = roipoly(I);                
        save(mask_file, 'BW');
        BW_rest(:,:,s) = BW;         
    end  
   end
   
    if(size(BW_stress,2)==size(BW_rest, 2))
        figure; imagescn(cat(4, BW_stress, BW_rest));
    else
        figure; imagescn(BW_stress);
        figure; imagescn(BW_rest);
    end
    clear BW_stress BW_rest
        
    pause;
    closeall;
end
