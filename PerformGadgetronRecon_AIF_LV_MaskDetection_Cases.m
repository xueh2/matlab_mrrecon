
function PerformGadgetronRecon_AIF_LV_MaskDetection_Cases(perf_cases, resDir, fig_Dir)
% PerformGadgetronRecon_AIF_LV_MaskDetection_Cases(perf_cases, fig_Dir)

startN = 1
endN = size(perf_cases, 1)

upslope_thres = 0.65
upslope_high_thres = 0.75 
auc_thres = 0.92
area_thres = 10
ro_boundary_ratio = 0.1
e1_boundary_ratio = 0.1

failed_cases = [];

for ii=startN:endN
    ii
    try
        disp([num2str(ii) ' of ' num2str(endN) ' ... ']);
        [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(perf_cases{ii,2});

        % stress
        stress_dir = fullfile(resDir, study_dates, perf_cases{ii,2})

        aif_moco_stress = analyze75read(fullfile(stress_dir, 'DebugOutput', 'aif_upsampled.hdr'));
        aif_pd_stress = analyze75read(fullfile(stress_dir, 'DebugOutput', 'aifPD_for_TwoEcho_T2StartCorrection_0.hdr'));
        aif_moco_stress_mask = analyze75read(fullfile(stress_dir, 'DebugOutput', 'aif_LV_mask.hdr'));

        [mask_stress, LV, peakTime_LV_stress] = PerformGadgetronRecon_AIF_LV_MaskDetection(aif_moco_stress, upslope_thres, auc_thres, area_thres, ro_boundary_ratio, e1_boundary_ratio);
       
        ind = find(mask_stress>0);
        aif_slc = aif_moco_stress(:,:, floor(peakTime_LV_stress));
        
        aif_slc_mask = aif_slc;
        aif_slc_mask(ind) = max(aif_slc(:))+1;
        
        ind = find(aif_moco_stress_mask>0);
        aif_slc_mask_old = aif_slc;
        aif_slc_mask_old(ind) = max(aif_slc(:))+1;
        
        h = figure('Name',perf_cases{ii,2},'NumberTitle','off'); imagescn(cat(3, aif_slc_mask, aif_slc, aif_slc_mask_old), [], [1 3], [10]);
        saveas(h, fullfile(fig_Dir, perf_cases{ii,2}), 'fig');
        ff = getframe(h);
        [X, map] = frame2im(ff);
        imwrite(X, fullfile(fig_Dir, [perf_cases{ii,2} '.bmp']), 'BMP');
        
        % rest
        if(size(perf_cases, 2)>=3)
            [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(perf_cases{ii, 3});
            rest_dir = fullfile(resDir, study_dates, perf_cases{ii,3})

            aif_moco_rest = analyze75read(fullfile(rest_dir, 'DebugOutput', 'aif_upsampled.hdr'));
            aif_pd_rest = analyze75read(fullfile(rest_dir, 'DebugOutput', 'aifPD_for_TwoEcho_T2StartCorrection_0.hdr'));
            aif_moco_rest_mask = analyze75read(fullfile(rest_dir, 'DebugOutput', 'aif_LV_mask.hdr'));

            [mask, LV, peakTime_LV] = PerformGadgetronRecon_AIF_LV_MaskDetection(aif_moco_rest, upslope_thres, auc_thres, area_thres, ro_boundary_ratio, e1_boundary_ratio);

            ind = find(mask>0);
            aif_slc = aif_moco_rest(:,:, floor(peakTime_LV));

            aif_slc_mask = aif_slc;
            aif_slc_mask(ind) = max(aif_slc(:))+1;

            ind = find(aif_moco_rest_mask>0);
            aif_slc_mask_old = aif_slc;
            aif_slc_mask_old(ind) = max(aif_slc(:))+1;

            h = figure('Name',perf_cases{ii,3},'NumberTitle','off'); imagescn(cat(3, aif_slc_mask, aif_slc, aif_slc_mask_old), [], [1 3], [10]);
            saveas(h, fullfile(fig_Dir, perf_cases{ii,3}), 'fig');
            ff = getframe(h);
            [X, map] = frame2im(ff);
            imwrite(X, fullfile(fig_Dir, [perf_cases{ii,3} '.bmp']), 'BMP');
        end
        
%         pause
        closeall
    catch
    end
end
