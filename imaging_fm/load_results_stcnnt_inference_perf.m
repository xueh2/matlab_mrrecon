function [input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, plot_flag)
% [input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, plot_flag)

input = [];
if exist(fullfile(res_dir, 'input.npy'))
    input = double(py.numpy.load(fullfile(res_dir, 'input.npy')));
else
    if exist(fullfile(res_dir, 'input_real.npy'))
        input = readNPY(fullfile(res_dir, 'input_real')) + i*readNPY(fullfile(res_dir, 'input_imag'));
    end
end

res = [];
if exist(fullfile(res_dir, 'output.npy'))
    res = double(py.numpy.load(fullfile(res_dir, 'output.npy')));    
else
    if exist(fullfile(res_dir, 'output_real.npy'))
        res = readNPY(fullfile(res_dir, 'output_real')) + j*readNPY(fullfile(res_dir, 'output_imag'));
    end
end

gmap = [];
if exist(fullfile(res_dir, 'gmap.npy'))
    gmap = readNPY(fullfile(res_dir, 'gmap.npy'));
end

if plot_flag
    if ndims(res) == 4
        SLC = size(res, 4);
        if SLC < 3
            h = figure('Name', res_dir); imagescn(abs(cat(4, input, res, input-res)), [0 4*mean(abs(res(:)))], [SLC 3], [8], 3);
        else
            h = figure('Name', res_dir); imagescn(abs(cat(4, input, res, input-res)), [0 4*mean(abs(res(:)))], [3 SLC], [8], 3);
        end
    end
    if ndims(res) == 5
        SLC = size(res, 5);
        h = figure('Name', res_dir); imagescn(abs(cat(4, input, res, input-res)), [0 4*mean(abs(res(:)))], [2*3 SLC], [8], 3);
    end
    if ndims(res) == 3
        h = figure('Name', res_dir); imagescn(abs(cat(4, input, res, input-res)), [], [1 3], [8], 3);
    end
end

size(input)
size(res)
