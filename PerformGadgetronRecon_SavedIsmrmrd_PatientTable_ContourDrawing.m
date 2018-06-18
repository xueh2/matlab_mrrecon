
function PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ContourDrawing(PerfTable, resDir, contourDir, SelectedType, splenic_column, selected_column, stress_column, rest_column, report_column)
% PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ContourDrawing(PerfTable, resDir, contourDir)
% PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ContourDrawing(PerfTable, 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area_Aif_recorded', 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area_ROI')
% PerformGadgetronRecon_SavedIsmrmrd_PatientTable_ContourDrawing(PerfTable, 'D:\data\ut\NewData\PaperResults\BARTS_Area_Aif_recorded', 'D:\data\ut\NewData\PaperResults\BARTS_Area_ROI', 20, 24)

if(nargin<4)
    SelectedType = [];
end

if(exist(contourDir)~=7)
    mkdir(contourDir);
end

if(isnumeric(PerfTable{1}))
    num = size(PerfTable, 1);
    start = 0;
else
    num = size(PerfTable, 1)-1;
    start = 1;
end

for n=start:num
    % disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, report_column} ' - ' num2str(PerfTable{n+1, rest_column+1}) ' - ' num2str(PerfTable{n+1, 17}) ' - ' num2str(PerfTable{n+1, 20})]); 
    disp(['=============================================================================================================================================']);  
    if(~isempty(report_column) & report_column>0)
        disp([num2str(n) ' out of ' num2str(num) ' - ' PerfTable{n+1, report_column} ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column} ' - ' PerfTable{n+1, splenic_column} ' - ' PerfTable{n+1, selected_column}]); 
    else
        disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]); 
    end
    PerfTable{n+1, :}
    
    stressCase = PerfTable{n+1, stress_column};
    restCase = PerfTable{n+1, rest_column};
    
    process = 0;
    if(isempty(SelectedType)~=1)
        
        selected = PerfTable{n+1, selected_column};

        for ii=1:numel(SelectedType)
            if( ~isempty(strfind(selected, SelectedType{ii})) )
                process = 1;
                disp(['Selected Type : ' selected]); 
                break;
            end
        end
    else
        process = 1;        
    end
    
%     if( (strcmp(splenic_cut_off, 'Yes') ~= 1) & (isempty( strfind(splenic_cut_off, 'yes') ) == 1) )
%         continue;
%     end
    
    if(~process)
        continue;
    end
        
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(stressCase);
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(restCase);
    
    roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' patientID '_' studyID '_' study_dates '_ROI'])
    
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
        
        resCaseDir = [scannerID '_' patientID '_' studyID '_' study_dates];    

        stressDir = fullfile(resDir, study_dates, stressCase);
        
        ind = find(stressDir=='\');
        ind2 = strfind(stressDir, scannerID);
        data_name = stressDir(ind(end)+1:ind2(1)-2);
        figDir = fullfile(resDir, study_dates, [data_name '_' resCaseDir '_' study_time_stress '_' study_time_rest '_Figure']);
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
    
    % AIF mask

%     restDir = fullfile(resDir, study_dates, restCase)
%     stressDir = fullfile(resDir, study_dates, stressCase)
% 
%     r1 = analyze75read(fullfile(restDir, 'DebugOutput', 'input_aif_moco_0.hdr'));
%     r2 = analyze75read(fullfile(restDir, 'DebugOutput', 'input_aif_moco_1.hdr'));
%     aif_moco_rest = cat(4, r1, r2);    
%     
%     r1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'input_aif_moco_0.hdr'));
%     r2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'input_aif_moco_1.hdr'));
%     aif_moco_stress = cat(4, r1, r2);    
%    
%     aif_rest_mask = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_LV_mask_for_TwoEcho_T2StartCorrection_0.hdr'));
%     aif_rest_mask_final = analyze75read(fullfile(restDir, 'DebugOutput', 'AifLVMask_after_Picking.hdr'));
% 
%     figure; imagescn(aif_moco_rest, [], [], [], 3);
%     figure; imagescn(aif_moco_stress, [], [], [], 3);
%         
%     AIF_stress_mask = load(fullfile(roiDir, 'AIF_stress_mask.mat'));
%     AIF_rest_mask = load(fullfile(roiDir, 'AIF_rest_mask.mat'));
%     
%     rep = size(r1, 3);
%     
%     BW=zeros(size(r1(:,:,1)));
%     BW=roipoly(r1(:,:,1),AIF_stress_mask.ROI_info_table(1).ROI_x_coordinates, AIF_stress_mask.ROI_info_table(1).ROI_y_coordinates);
%     BW_stress_AIF = flipdim(BW,2);
%     
%     BW=roipoly(r1(:,:,1),AIF_rest_mask.ROI_info_table(1).ROI_x_coordinates, AIF_rest_mask(1).ROI_info_table(1).ROI_y_coordinates);
%     BW_rest_AIF = flipdim(BW,2);
%     
%     BWs = repmat(BW_stress_AIF, [1 1 rep]);
%     BWr = repmat(BW_rest_AIF, [1 1 rep]);
%     
%     ind = find(BWs>0);
%     ind2 = find(BWr>0);
%     
%     r1(ind) = 1024;
%     r2(ind2) = 1024;
%     aif_moco_stress = cat(4, r1, r2); 
%     figure; imagescn(aif_moco_stress, [], [], 10, 3);
    
%     header = CreateGtImageHeader(BW_stress_AIF);
%     Matlab_gt_write_analyze(single(BW_stress_AIF), header, fullfile(roiDir, 'AIF_stress_mask'));
%     Matlab_gt_write_analyze(single(BW_rest_AIF), header, fullfile(roiDir, 'AIF_rest_mask'));
    
    pause;
    closeall;
end
