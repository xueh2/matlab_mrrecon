
function PerformGadgetronRecon_Perf_SpatialTemporalFiltering(resDir, perf_cases, thres_temporal, thres_spatial, wave_name, level)
% PerformGadgetronRecon_Perf_SpatialTemporalFiltering(resDir, perf_cases, thres_temporal, thres_spatial, wave_name, level)

num = size(perf_cases, 1);

suffix = ['temporal_' num2str(thres_temporal) '_spatial_' num2str(thres_spatial) '_wave_' wave_name '_level_' num2str(level)];
ind = find(suffix=='.');
suffix(ind) = 'p';

for n=1:num
    
    name = perf_cases{n, 1};
    disp(['Process ' num2str(n) ' out of ' num2str(num) ' - ' name]);
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
    
    res_dir_case = fullfile(resDir, study_dates, name);
    
    if(~isFileExist(res_dir_case))
        disp(['Case is not processed : ' name]);
        continue;
    end    
    
    debug_dir = fullfile(res_dir_case, 'DebugOutput');
    
    [ori, filtered] = PerformGadgetronRecon_Perf_Debug_SpatialTemporalFiltering(debug_dir, thres_temporal, thres_spatial, wave_name, level);
    
%     SLC = 1;
%     while(isFileExist(fullfile(debug_dir, ['perf_moco_' num2str(SLC-1) '.hdr'])))
%         SLC = SLC + 1;
%     end
%     
%     for slc=1:SLC-1
%         data = analyze75read( fullfile(debug_dir, ['perf_moco_' num2str(slc-1) '.hdr']) );
%         res = PerformSpatioTemporalFiltering(complex(double(data)), thres_temporal, thres_spatial, wave_name, level);        
%         header = CreateGtImageHeader(data, [1 1 1 1]);
%         Matlab_gt_write_analyze(single(abs(res)), header, fullfile(debug_dir, ['perf_moco_' num2str(slc-1) '_' suffix '.hdr']));
%     end    
end
