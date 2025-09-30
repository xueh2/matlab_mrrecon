function [ssims, psnrs, noise_sigmas, gmap, clean_im, input, outputs] = compare_models_metric(case_dir, models_dir, plot_flag, only_metrics)
% [ssims, psnrs, noise_sigmas, gmap, clean_im, outputs] =
% compare_models_metric('/data/raw_data/retro_cine_3T/raw_selected/00047850/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071167_591071176_821_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE/model/res/slc_0/', {'hrnet_28m', 'hrnet_50m', 'hrnet_100m'}, 1);

N = numel(models_dir);

noise_sigmas = readNPY(fullfile(case_dir, 'noise_sigmas.npy'));
gmap = readNPY(fullfile(case_dir, 'gmap.npy'));

clean_im = [];
input = [];
outputs = [];

if ~only_metrics
    res_R2 = fullfile(case_dir, '../../res_R2');
    if exist(res_R2)
        disp(['Loading from res R2: ' res_R2])
        clean_im = complex(readNPY(fullfile(res_R2, 'output_real')), readNPY(fullfile(res_R2, 'output_imag')));
        clean_im = squeeze(clean_im);
    else
        clean_im = complex(readNPY(fullfile(case_dir, '../../input_real')), readNPY(fullfile(case_dir, '../../input_imag'))); 
        clean_im = squeeze(clean_im);
    end
    input = complex(readNPY(fullfile(case_dir, 'input_real')), readNPY(fullfile(case_dir, 'input_imag')));
end

ind = strfind(case_dir, '/slc_');
slc = str2num(case_dir(ind+5:end));
disp([case_dir ' - slc is ' num2str(slc)])

if ~only_metrics
    if size(clean_im, 4) > 1
        clean_im = clean_im(:,:,:,slc+1);
    end
    output = complex(readNPY(fullfile(case_dir, models_dir{1}, 'output_real')), readNPY(fullfile(case_dir, models_dir{1}, 'output_imag')));
    
    RO = size(output, 1);
    E1 = size(output, 2);
    PHS = size(output, 3);
    num_snr_levels = size(output, 4);

    outputs = zeros(RO, E1, PHS, num_snr_levels, N);
    a = squeeze(clean_im);

    ssims = zeros(N, num_snr_levels);
    psnrs = zeros(N, num_snr_levels);
else
    outputs = [];
    ssims = [];
    psnrs = [];
end

for n=1:N

    if ~exist(fullfile(case_dir, models_dir{n}, 'output_real.npy')) & ~exist(fullfile(case_dir, models_dir{n}, 'output.npy'))
        continue
    end
    disp(models_dir{n})

    if only_metrics & exist(fullfile(case_dir, models_dir{n}, 'metrics.mat'))
        rr = load(fullfile(case_dir, models_dir{n}, 'metrics.mat'));
        if isempty(ssims)
            num_snr_levels = size(rr.ssim_case, 2);
            ssims = zeros(N, num_snr_levels);
            psnrs = zeros(N, num_snr_levels);
        end
        ssims(n, :) = rr.ssim_case;
        psnrs(n, :) = rr.psnr_case;
    else
        if exist(fullfile(case_dir, models_dir{n}, 'output_real.npy'))
            output = complex(readNPY(fullfile(case_dir, models_dir{n}, 'output_real')), readNPY(fullfile(case_dir, models_dir{n}, 'output_imag')));
        else
            output = readNPY(fullfile(case_dir, models_dir{n}, 'output'));
        end
    
        for nn=1:num_snr_levels
            
            b = double(output(:,:,:,nn) * noise_sigmas(nn));
            a_psnr = compute_psnr(abs(a), abs(b));
            a_ssim = multissim3(abs(a), abs(b));
    
            ssims(n, nn) = a_ssim;
            psnrs(n, nn) = a_psnr;
    
            outputs(:,:,:,nn, n) = b;
        end
    
        ssim_case = ssims(n, :);
        psnr_case = psnrs(n, :);
    
        save(fullfile(case_dir, models_dir{n}, 'metrics.mat'), 'ssim_case', 'psnr_case');
    end
end

if ~only_metrics
    for nn=1:num_snr_levels
        input(:,:,:,nn) = input(:,:,:,nn) * noise_sigmas(nn);
    end
end

if plot_flag
    clean_ims = repmat(squeeze(clean_im), [1 1 1 num_snr_levels]);
    h = figure('Name', case_dir);
    S = abs(cat(5, clean_ims, input, outputs));
    imagescn(reshape(S, [RO, E1, PHS, (N+2)*num_snr_levels]), [0 8*abs(mean(clean_im(:)))], [N+2, num_snr_levels], [16], 3);
end