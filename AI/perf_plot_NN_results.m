
function [h_Gd, h_fmap, h_fSD] = perf_plot_NN_results(trainingDataDir, endo, epi, endo_epi_rv_rvi)

data.fmap_resized_training = readNPY(fullfile(trainingDataDir, 'fmap_resized_training_norm.npy'));
data.Gd_resized_training = readNPY(fullfile(trainingDataDir, 'Gd_resized_training_norm.npy'));
data.fSD_resized_training = readNPY(fullfile(trainingDataDir, 'fSD_resized_training_norm.npy'));

SLC = size(data.Gd_resized_training, 4);

if(SLC==3)
    figure; imagescn(cat(3, endo.seg0, endo.seg1, endo.seg2, epi.seg0, epi.seg1, epi.seg2), [], [2 3], 12);

    h = figure; 
    imagescn(cat(3, endo.prob0, endo.prob1, endo.prob2, ... 
        epi.prob0, epi.prob1, epi.prob2, ... 
        endo_epi_rv_rvi.prob0_rv, endo_epi_rv_rvi.prob1_rv, endo_epi_rv_rvi.prob2_rv, ...
        endo_epi_rv_rvi.prob0_rvi, endo_epi_rv_rvi.prob1_rvi, endo_epi_rv_rvi.prob2_rvi), [], [4 3], 12);PerfColorMap;

    h_Gd = perf_plot_contours_on_images(data.Gd_resized_training, [{endo.C0}, {endo.C1}, {endo.C2}], [{epi.C0}, {epi.C1}, {epi.C2}], [{endo_epi_rv_rvi.rv0}, {endo_epi_rv_rvi.rv1}, {endo_epi_rv_rvi.rv2}], endo_epi_rv_rvi.rvi, [0 1.5]);
    h_fmap = perf_plot_contours_on_images(data.fmap_resized_training, [{endo.C0}, {endo.C1}, {endo.C2}], [{epi.C0}, {epi.C1}, {epi.C2}], [{endo_epi_rv_rvi.rv0}, {endo_epi_rv_rvi.rv1}, {endo_epi_rv_rvi.rv2}], endo_epi_rv_rvi.rvi, [0 8]);
    h_fSD = perf_plot_contours_on_images(data.fSD_resized_training, [{endo.C0}, {endo.C1}, {endo.C2}], [{epi.C0}, {epi.C1}, {epi.C2}], [{endo_epi_rv_rvi.rv0}, {endo_epi_rv_rvi.rv1}, {endo_epi_rv_rvi.rv2}], endo_epi_rv_rvi.rvi, [0 1]);
else
    figure; imagescn(cat(3, endo.seg0, endo.seg1, endo.seg2, endo.seg3, epi.seg0, epi.seg1, epi.seg2, epi.seg3), [], [2 4], 12);

    h = figure; 
    imagescn(cat(3, endo.prob0, endo.prob1, endo.prob2, endo.prob3, ... 
        epi.prob0, epi.prob1, epi.prob2, epi.prob3, ... 
        endo_epi_rv_rvi.prob0_rv, endo_epi_rv_rvi.prob1_rv, endo_epi_rv_rvi.prob2_rv, endo_epi_rv_rvi.prob3_rv, ...
        endo_epi_rv_rvi.prob0_rvi, endo_epi_rv_rvi.prob1_rvi, endo_epi_rv_rvi.prob2_rvi, endo_epi_rv_rvi.prob3_rvi), [], [4 4], 12);PerfColorMap;

    h_Gd = perf_plot_contours_on_images(data.Gd_resized_training, [{endo.C0}, {endo.C1}, {endo.C2}, {endo.C3}], [{epi.C0}, {epi.C1}, {epi.C2}, {epi.C3}], [{endo_epi_rv_rvi.rv0}, {endo_epi_rv_rvi.rv1}, {endo_epi_rv_rvi.rv2}, {endo_epi_rv_rvi.rv3}], endo_epi_rv_rvi.rvi, [0 1.5]);
    h_fmap = perf_plot_contours_on_images(data.fmap_resized_training, [{endo.C0}, {endo.C1}, {endo.C2}, {endo.C3}], [{epi.C0}, {epi.C1}, {epi.C2}, {epi.C3}], [{endo_epi_rv_rvi.rv0}, {endo_epi_rv_rvi.rv1}, {endo_epi_rv_rvi.rv2}, {endo_epi_rv_rvi.rv3}], endo_epi_rv_rvi.rvi, [0 8]);
    h_fSD = perf_plot_contours_on_images(data.fSD_resized_training, [{endo.C0}, {endo.C1}, {endo.C2}, {endo.C3}], [{epi.C0}, {epi.C1}, {epi.C2}, {epi.C3}], [{endo_epi_rv_rvi.rv0}, {endo_epi_rv_rvi.rv1}, {endo_epi_rv_rvi.rv2}, {endo_epi_rv_rvi.rv3}], endo_epi_rv_rvi.rvi, [0 1]);
end

