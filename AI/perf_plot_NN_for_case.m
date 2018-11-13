
function [data, endo, epi, endo_epi_rv_rvi] = perf_plot_NN_for_case(resDir, caseName, only_view)
% [data, endo, epi, endo_epi_rv_rvi] = perf_plot_NN_for_case(resDir, caseName, only_view)

[configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(caseName);
trainingDataDir = fullfile(resDir, study_date, caseName, 'Training_Data');

data_file = fullfile(trainingDataDir, ['training_data.mat'])
try
    data = load(data_file);
catch
    disp(['cannot load ' data_file])
    failed_cases = [failed_cases; {caseName}];
end

endo_model = 'perf_endo_network_2018-11-05.pbt'
epi_model = 'perf_epi_network_2018-11-05.pbt'
endo_epi_rv_rvi_model = 'perf_endo_epi_rv_rvi_network_2018-11-04.pbt'

if(only_view)
    [endo, epi, endo_epi_rv_rvi] = perf_apply_NN_on_images(trainingDataDir, [], [], []);    
else
    [endo, epi, endo_epi_rv_rvi] = perf_apply_NN_on_images(trainingDataDir, endo_model, epi_model, endo_epi_rv_rvi_model);
end

figure; imagescn(data.fmap_resized_training, [0 8], [1 3], 16); PerfColorMap;

% ---------------------------------------------------------------------------------------------------------
% perf_add_endo_epi_contours(data.Gd_resized_training, trainingDataDir)

[h_Gd, h_fmap, h_fSD] = perf_plot_NN_results(trainingDataDir, endo, epi, endo_epi_rv_rvi);
