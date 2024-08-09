function [input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, plot_flag)
% [input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, plot_flag)

input = readNPY(fullfile(res_dir, 'input_real')) + i*readNPY(fullfile(res_dir, 'input_imag'));
res = readNPY(fullfile(res_dir, 'output_real')) + j*readNPY(fullfile(res_dir, 'output_imag'));
gmap = readNPY(fullfile(res_dir, 'gmap'));

if plot_flag
    SLC = size(input, 4);
    figure; imagescn(cat(4, input, res), [], [2 SLC], [16], 3);
end