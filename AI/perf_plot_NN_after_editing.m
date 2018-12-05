
function [h_Gd, h_fmap, h_fSD] = perf_plot_NN_after_editing(trainingDataDir, Seg)
% [h_Gd, h_fmap, h_fSD] = perf_plot_NN_after_editing(trainingDataDir, Seg)

endo_epi_rv_rvi_resized_training_mask = cat(4, Seg(1).endo_epi_rv_rvi_resized_training_mask, Seg(2).endo_epi_rv_rvi_resized_training_mask, Seg(3).endo_epi_rv_rvi_resized_training_mask);
figure; imagescn(endo_epi_rv_rvi_resized_training_mask, [], [1 3], 12);

endo.C0 = Seg(1).endo_resized_training(:, 2:-1:1);
endo.C1 = Seg(2).endo_resized_training(:, 2:-1:1);
endo.C2 = Seg(3).endo_resized_training(:, 2:-1:1);

epi.C0 = Seg(1).epi_resized_training(:, 2:-1:1);
epi.C1 = Seg(2).epi_resized_training(:, 2:-1:1);
epi.C2 = Seg(3).epi_resized_training(:, 2:-1:1);

endo_epi_rv_rvi.rv0 = Seg(1).rv_resized_training(:, 2:-1:1);
endo_epi_rv_rvi.rv1 = Seg(2).rv_resized_training(:, 2:-1:1);
endo_epi_rv_rvi.rv2 = Seg(3).rv_resized_training(:, 2:-1:1);

endo_epi_rv_rvi.rvi = zeros(3,2);
endo_epi_rv_rvi.rvi(1,:) = Seg(1).rvi_resized_training(2:-1:1);
endo_epi_rv_rvi.rvi(2,:) = Seg(2).rvi_resized_training(2:-1:1);
endo_epi_rv_rvi.rvi(3,:) = Seg(3).rvi_resized_training(2:-1:1);

data.fmap_resized_training = readNPY(fullfile(trainingDataDir, 'fmap_resized_training.npy'));
data.Gd_resized_training = readNPY(fullfile(trainingDataDir, 'Gd_resized_training.npy'));
data.fSD_resized_training = readNPY(fullfile(trainingDataDir, 'fSD_resized_training.npy'));

h_Gd = perf_plot_contours_on_images(data.Gd_resized_training, [{endo.C0}, {endo.C1}, {endo.C2}], [{epi.C0}, {epi.C1}, {epi.C2}], [{endo_epi_rv_rvi.rv0}, {endo_epi_rv_rvi.rv1}, {endo_epi_rv_rvi.rv2}], endo_epi_rv_rvi.rvi, [0 1.5]);
h_fmap = perf_plot_contours_on_images(data.fmap_resized_training, [{endo.C0}, {endo.C1}, {endo.C2}], [{epi.C0}, {epi.C1}, {epi.C2}], [{endo_epi_rv_rvi.rv0}, {endo_epi_rv_rvi.rv1}, {endo_epi_rv_rvi.rv2}], endo_epi_rv_rvi.rvi, [0 8]);
h_fSD = perf_plot_contours_on_images(data.fSD_resized_training, [{endo.C0}, {endo.C1}, {endo.C2}], [{epi.C0}, {epi.C1}, {epi.C2}], [{endo_epi_rv_rvi.rv0}, {endo_epi_rv_rvi.rv1}, {endo_epi_rv_rvi.rv2}], endo_epi_rv_rvi.rvi, [0 1]);
