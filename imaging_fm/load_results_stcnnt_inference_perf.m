function [input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, plot_flag)
% [input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, plot_flag)

input = [];
if exist(fullfile(res_dir, 'input_real.npy'))
    input = readNPY(fullfile(res_dir, 'input_real')) + i*readNPY(fullfile(res_dir, 'input_imag'));
end
if exist(fullfile(res_dir, 'input.npy'))
    input = double(py.numpy.load(fullfile(res_dir, 'input.npy')));
end

if exist(fullfile(res_dir, 'output_real.npy'))
    res = readNPY(fullfile(res_dir, 'output_real')) + j*readNPY(fullfile(res_dir, 'output_imag'));
else
    res = double(py.numpy.load(fullfile(res_dir, 'output.npy')));
end

gmap = [];
if exist(fullfile(res_dir, 'gmap.npy'))
    gmap = readNPY(fullfile(res_dir, 'gmap'));
end

if plot_flag
    SLC = size(res, 4);
    if SLC < 3
        figure; imagescn(abs(cat(4, input, res, input-res)), [0 4*mean(abs(res(:)))], [SLC 3], [8], 3);
    else
        figure; imagescn(abs(cat(4, input, res, input-res)), [0 4*mean(abs(res(:)))], [3 SLC], [8], 3);
    end
end

size(input)
size(res)
