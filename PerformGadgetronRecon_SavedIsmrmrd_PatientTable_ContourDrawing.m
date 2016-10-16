
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

num = size(PerfTable, 1)-1;
for n=1:num
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column} ' - ' PerfTable{n+1, splenic_column} ' - ' PerfTable{n+1, selected_column}]); 
    % disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, report_column} ' - ' num2str(PerfTable{n+1, rest_column+1}) ' - ' num2str(PerfTable{n+1, 17}) ' - ' num2str(PerfTable{n+1, 20})]); 
    disp(['==================================================================']);  
    
    splenic_cut_off = PerfTable{n+1, splenic_column};
    selected = PerfTable{n+1, selected_column};
    stressCase = PerfTable{n+1, stress_column};
    restCase = PerfTable{n+1, rest_column};
    
    process = 0;
    if(isempty(SelectedType)~=1)
        for ii=1:numel(SelectedType)
            if( ~isempty(strfind(selected, SelectedType{ii})) )
                process = 1;
                break;
            end
        end
    else
        process = 1;
    end
    
    if( (strcmp(splenic_cut_off, 'Yes') ~= 1) & (isempty( strfind(splenic_cut_off, 'yes') ) == 1) )
        continue;
    end
    
    if(~process)
        continue;
    end
    
    disp(['Selected Type : ' selected]); 
    
    [h_flow_stress, h_flow_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir,  stressCase, restCase, 1);        
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(stressCase);
    
    roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' patientID '_' studyID '_' study_dates '_ROI']);
    
    if(isFileExist(roiDir))
        dir(roiDir);        
    else    
        mkdir(roiDir)
    end
    cd(roiDir)
    
    for s=1:numel(h_flow_stress)
        [I, xlim, ylim, clim, position]=getimagedata(h_flow_stress(s));        
        BW = roipoly(I);        
        mask_file = fullfile(roiDir, ['perf_mask_stress_' num2str(s-1) '.mat']);
        save(mask_file, 'BW');
    end
    
    for s=1:numel(h_flow_rest)
        [I, xlim, ylim, clim, position]=getimagedata(h_flow_rest(s));        
        BW = roipoly(I);        
        mask_file = fullfile(roiDir, ['perf_mask_rest_' num2str(s-1) '.mat']);
        save(mask_file, 'BW');
    end
    
    pause;
    closeall;
end
