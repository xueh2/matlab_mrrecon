
function PerformGadgetronRecon_SavedIsmrmrd_Perfusion_Matlab(PerfTable, dataDir, resDir, roiDir, reComputed, with_hct, splenic_column, selected_column, stress_column, rest_column, hct_column, hct_preset, reComputed_OnlyGlobal, reComputed_withoutR2Star, case_list, res_string, res_string_Q_e)
% PerformGadgetronRecon_SavedIsmrmrd_Perfusion_Matlab(PerfTable, dataDir, resDir, roiDir, reComputed, with_hct, splenic_column, selected_column, stress_column, rest_column, hct_column, hct_preset, reComputed_OnlyGlobal, reComputed_withoutR2Star)
% PerformGadgetronRecon_SavedIsmrmrd_Perfusion_Matlab(PerfTable, resDir, reComputed, splenic_column, selected_column, stress_column, rest_column, hct_column)
% PerformGadgetronRecon_SavedIsmrmrd_Perfusion_Matlab(PerfTable, 'I:\KAROLINSKA', 'I:\ReconResults\KAROLINSKA')

if(nargin<11)
    hct_preset = -1;
end

if(nargin<12)
    reComputed_OnlyGlobal = 0;
end

if(nargin<13)
    reComputed_withoutR2Star = 0;
end

num = size(PerfTable, 1)-1;

if(nargin<15 | isempty(case_list))
    case_list = 1:num;
end

if(nargin<16)
    res_string = 'Linear';
end

if(nargin<17)
    res_string_Q_e = res_string;
end

for tt=1:numel(case_list)
    
    n = case_list(tt);
    
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column} ' - ' PerfTable{n+1, splenic_column} ' - ' PerfTable{n+1, selected_column}]); 
    disp(['==================================================================']);  
    
    splenic_cut_off = PerfTable{n+1, splenic_column};
    selected = PerfTable{n+1, selected_column};
    stressCase = PerfTable{n+1, stress_column};
    restCase = PerfTable{n+1, rest_column};
    
    if(hct_preset>=0 & hct_preset<=1)
        HCT = hct_preset;
    else
        HCT = PerfTable{n+1, hct_column};
        if(~with_hct)
            HCT=0;    
        end
    end
           
    if( (strcmp(splenic_cut_off, 'Yes') ~= 1) & (isempty( strfind(splenic_cut_off, 'yes') ) == 1) )
        continue;
    end
          
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(stressCase);
        
    if(isnan(HCT))
        HCT = 0;
    else
        if(HCT>1)
            HCT = HCT/100;
        end
    end
    
    roiDir_case = fullfile(roiDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' patientID '_' studyID '_' study_dates '_ROI'])
    
    % stress
    try
    PerformGadgetronRecon_Matlab_FlowMapping_Linear(dataDir, stressCase, resDir, roiDir_case, 'stress', HCT, reComputed, reComputed_OnlyGlobal, reComputed_withoutR2Star, res_string, res_string_Q_e);
    catch
    end
    
    %rest
    try
    PerformGadgetronRecon_Matlab_FlowMapping_Linear(dataDir, restCase, resDir, roiDir_case, 'rest', HCT, reComputed, reComputed_OnlyGlobal, reComputed_withoutR2Star, res_string, res_string_Q_e);
    catch
    end
%     pause;
    closeall;
end
