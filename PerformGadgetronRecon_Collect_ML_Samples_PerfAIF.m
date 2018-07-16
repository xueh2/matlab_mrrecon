
function [error_lists, weak_aif_lists] = PerformGadgetronRecon_Collect_ML_Samples_PerfAIF(perf_cases, rest_cases, dataDir, resDir, host, configNamePreset, sample_dir, run_case)
% PerformGadgetronRecon_Collect_ML_Samples_PerfAIF(perf_cases, rest_cases, dataDir, resDir, host, configNamePreset, sample_dir)

case_lists = []
error_lists = []
weak_aif_lists = []

startN = 1
endN = size(perf_cases, 1)
for ii=startN:endN   
    case_lists = [case_lists; perf_cases(ii, 2); perf_cases(ii, 3) ];
end

endN = size(rest_cases, 1)
try
    for ii=startN:endN   
        case_lists = [case_lists; rest_cases(ii, 2);];
    end
catch
    for ii=startN:endN   
        case_lists = [case_lists; rest_cases(ii, 1);];
    end
end

disp(['Total number of cases : ' num2str(numel(case_lists))])

for ii=startN:numel(case_lists)   

    disp([num2str(ii) '/' num2str(numel(case_lists)) ' - ' case_lists{ii}])
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(case_lists{ii});
    case_dir = fullfile(resDir, study_dates, case_lists{ii});        
    h5Name = fullfile(dataDir, study_dates, [case_lists{ii} '.h5']);
    
    if (isFileExist(fullfile(sample_dir, case_lists{ii}, 'aif.mat')))
        continue;
    end
    
    if (~isempty(strfind(case_lists{ii}, 'Copy')))
        continue;
    end   
    
    try
        aif_signal = analyze75read(fullfile(case_dir, 'DebugOutput', 'aif_cin_echo0_all_signal_baseline_corrected.hdr'));
        
        if(max(aif_signal)<1.2)
            disp('max(aif_cin_Gd_baseline_corrected)<1.2')
            weak_aif_lists = [weak_aif_lists; case_lists(ii)];
%             continue;
        end
        
        save_case(case_dir, sample_dir, case_lists{ii});
        
%         [perf, ori, moco, moco_norm, PD, input_for_filter, filtered, aif_acq_time, perf_acq_time, dst_acq_time, ... 
%         aif_im, aif_moco, aif_cin, aif_cin_Gd, aif_cin_Gd_without_R2Star, aif_cin_Gd_baseline_corrected, ... 
%         aif_cin_all_echo0_signal, aif_cin_all_echo1_signal, aif_cin_all_echo0_signal_after_R2StarCorrection, ...
%         aif_cin_all_echo0_OverPD_after_R2StarCorrection, aif_cin_all_R2Star,  aif_cin_all_R2Star_SLEP, ... 
%         aif_PD, aif_mask, aif_mask_final, aif_LV_mask_plot, aif, aif_baseline_corrected, aif_plots, ... 
%         flow, Ki, PS, Vp, Visf, E, SDMap, Delay, ...
%         BTEX_Flow_all, BTEX_PS_all, BTEX_Vp_all, BTEX_Visf_all, BTEX_cost_all, ... 
%         BTEX_flow_SD_all, BTEX_PS_SD_all, BTEX_Visf_SD_all, BTEX_Vp_SD_all, BTEX_cov_all, ...
%         flow_SD, PS_SD, Vp_SD, Visf_SD, BTEX_cov, ... 
%         CC_F_PS, CC_F_Vp, CC_F_Visf, CC_PS_Vp, CC_PS_Visf, CC_Vp_Visf, ... 
%         BTEX_Tc_all, Fermi_Delay] = read_in_GT_Perf_DebugOutput_results(case_dir);
% 
%         aif_moco_echo1 = aif_moco(:,:,:,1);
%     
%         mkdir(fullfile(sample_dir, case_lists{ii}));        
%         cd(fullfile(sample_dir, case_lists{ii}))
% 
%         save perf aif_moco_echo1 ori moco moco_norm PD input_for_filter filtered aif_acq_time perf_acq_time dst_acq_time aif_im aif_moco aif_cin aif_cin_Gd aif_cin_Gd_without_R2Star aif_cin_Gd_baseline_corrected aif_cin_all_echo0_signal aif_cin_all_echo1_signal aif_cin_all_echo0_signal_after_R2StarCorrection aif_cin_all_echo0_OverPD_after_R2StarCorrection aif_cin_all_R2Star  aif_cin_all_R2Star_SLEP aif_PD aif_mask aif_mask_final aif_LV_mask_plot aif aif_baseline_corrected aif_plots flow Ki PS Vp Visf E SDMap Delay BTEX_Flow_all BTEX_PS_all BTEX_Vp_all BTEX_Visf_all BTEX_cost_all BTEX_flow_SD_all BTEX_PS_SD_all BTEX_Visf_SD_all BTEX_Vp_SD_all BTEX_cov_all flow_SD PS_SD Vp_SD Visf_SD BTEX_cov CC_F_PS CC_F_Vp CC_F_Visf CC_PS_Vp CC_PS_Visf CC_Vp_Visf BTEX_Tc_all Fermi_Delay
%         
%         save aif aif_moco_echo1 aif_acq_time aif_im aif_moco aif_cin aif_cin_Gd aif_cin_Gd_without_R2Star aif_cin_Gd_baseline_corrected aif_cin_all_echo0_signal aif_cin_all_echo1_signal aif_cin_all_echo0_signal_after_R2StarCorrection aif_cin_all_echo0_OverPD_after_R2StarCorrection aif_cin_all_R2Star  aif_cin_all_R2Star_SLEP aif_PD aif_mask aif_mask_final aif_LV_mask_plot aif aif_baseline_corrected aif_plots 
    catch
        
        if(~run_case)
            continue;
        end
        
        perf_xml = find_perf_xml(configNamePreset, case_lists{ii});
        checkProcessed = 0
        delete_old_res = 1
        startRemoteGT = 1
        PerformGadgetronRecon_SavedIsmrmrd_OneType_OneData(dataDir, case_lists(ii), host, resDir, checkProcessed, delete_old_res, startRemoteGT, {perf_xml});
        
        error_lists = [error_lists; case_lists(ii)];
        
        try
            save_case(case_dir, sample_dir, case_lists{ii});
        catch
            disp([case_lists{ii} ' should be removed ... '])
        end
    end
end
end

function save_case(case_dir, sample_dir, case_name)

    only_aif = 1;

    [perf, ori, moco, moco_norm, PD, input_for_filter, filtered, aif_acq_time, perf_acq_time, dst_acq_time, ... 
        aif_im, aif_moco, aif_cin, aif_cin_Gd, aif_cin_Gd_without_R2Star, aif_cin_Gd_baseline_corrected, ... 
        aif_cin_all_echo0_signal, aif_cin_all_echo1_signal, aif_cin_all_echo0_signal_after_R2StarCorrection, ...
        aif_cin_all_echo0_OverPD_after_R2StarCorrection, aif_cin_all_R2Star,  aif_cin_all_R2Star_SLEP, ... 
        aif_PD, aif_mask, aif_mask_final, aif_LV_mask_plot, aif, aif_baseline_corrected, aif_plots, ... 
        flow, Ki, PS, Vp, Visf, E, SDMap, Delay, ...
        BTEX_Flow_all, BTEX_PS_all, BTEX_Vp_all, BTEX_Visf_all, BTEX_cost_all, ... 
        BTEX_flow_SD_all, BTEX_PS_SD_all, BTEX_Visf_SD_all, BTEX_Vp_SD_all, BTEX_cov_all, ...
        flow_SD, PS_SD, Vp_SD, Visf_SD, BTEX_cov, ... 
        CC_F_PS, CC_F_Vp, CC_F_Visf, CC_PS_Vp, CC_PS_Visf, CC_Vp_Visf, ... 
        BTEX_Tc_all, Fermi_Delay] = read_in_GT_Perf_DebugOutput_results(case_dir, only_aif);

        aif_moco_echo1 = aif_moco(:,:,:,1);

        mkdir(fullfile(sample_dir, case_name));        
        cd(fullfile(sample_dir, case_name))

        % save perf aif_moco_echo1 ori moco moco_norm PD input_for_filter filtered aif_acq_time perf_acq_time dst_acq_time aif_im aif_moco aif_cin aif_cin_Gd aif_cin_Gd_without_R2Star aif_cin_Gd_baseline_corrected aif_cin_all_echo0_signal aif_cin_all_echo1_signal aif_cin_all_echo0_signal_after_R2StarCorrection aif_cin_all_echo0_OverPD_after_R2StarCorrection aif_cin_all_R2Star  aif_cin_all_R2Star_SLEP aif_PD aif_mask aif_mask_final aif_LV_mask_plot aif aif_baseline_corrected aif_plots flow Ki PS Vp Visf E SDMap Delay BTEX_Flow_all BTEX_PS_all BTEX_Vp_all BTEX_Visf_all BTEX_cost_all BTEX_flow_SD_all BTEX_PS_SD_all BTEX_Visf_SD_all BTEX_Vp_SD_all BTEX_cov_all flow_SD PS_SD Vp_SD Visf_SD BTEX_cov CC_F_PS CC_F_Vp CC_F_Visf CC_PS_Vp CC_PS_Visf CC_Vp_Visf BTEX_Tc_all Fermi_Delay

        save aif aif_moco_echo1 aif_acq_time aif_im aif_moco aif_cin aif_cin_Gd aif_cin_Gd_without_R2Star aif_cin_Gd_baseline_corrected aif_cin_all_echo0_signal aif_cin_all_echo1_signal aif_cin_all_echo0_signal_after_R2StarCorrection aif_cin_all_echo0_OverPD_after_R2StarCorrection aif_cin_all_R2Star  aif_cin_all_R2Star_SLEP aif_PD aif_mask aif_mask_final aif_LV_mask_plot aif aif_baseline_corrected aif_plots 
end
