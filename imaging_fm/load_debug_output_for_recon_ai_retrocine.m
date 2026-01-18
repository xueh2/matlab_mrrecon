function [kspace, aliasedIm, acs_src, acs_dst, coil_map, recon_kspace, recon_im, ker, kIm, unmixing, gmap] = load_debug_output_for_recon_ai_retrocine(debug_folder, accel, plot_flag)
% [kspace, ref_src, ref_dst, coil_map, recon_kspace, recon_im, ker, kIm, unmixing, gmap] = load_debug_output_for_recon_ai_retrocine(debug_folder, plot_flag)

kspace = readGTPlusExportData(fullfile(debug_folder, 'data_src_encoding_0'));
size(kspace)

RO = size(kspace, 1);
E1 = size(kspace, 2);
E2 = size(kspace, 3);
CHA = size(kspace, 4);
N = size(kspace, 5);
S = size(kspace, 6);
SLC = size(kspace, 7);
kspace = reshape(kspace, [RO, E1, E2, CHA, N, S, SLC]);

coil_map = readGTPlusExportData(fullfile(debug_folder, 'coil_map_encoding_0'));
size(coil_map)
coil_map = reshape(coil_map, [RO, E1, E2, CHA, 1, 1, SLC]);

recon_im = readGTPlusExportData(fullfile(debug_folder, 'recon_res_encoding_0_1'));
size(recon_im)
recon_im = reshape(recon_im, [RO, E1, E2, 1, N, 1, SLC]);

aliasedIm = readGTPlusExportData(fullfile(debug_folder, 'aliasedIm_encoding_0'));
size(aliasedIm)
aliasedIm = reshape(aliasedIm, [RO, E1, E2, CHA, N, 1, SLC]);

unWarppedIm = aliasedIm;

for slc=1:SLC
    a_kIm = readGTPlusExportData(fullfile(debug_folder, ['compute_kernel_2d_kIm_' num2str(slc-1)]));
    size(a_kIm)
    
    for n=1:N
        a_aliasedIm = squeeze(aliasedIm(:, :, 1, :, n,:,slc));
        unWarppedIm(:,:,:,:,n,1,slc) = sum(repmat(a_aliasedIm, [1 1 1 CHA]) .* a_kIm, 3);
    
        % for cprime = 1:CHA
        %     unaliasedIm = zeros(RO, E1);
        %     %a_aliasedIm = aliasedIm(:, :, :, n);
        %     for c=1:CHA
        %         unaliasedIm = unaliasedIm + aliasedIm(:, :, 1, c, n) .* a_kIm(:,:,c,cprime);
        %     end
        %     unWarppedIm(:,:,cprime, n) = unaliasedIm;
        % end
    end
end

recon_kspace = fft2c(unWarppedIm);

acs_src = [];
acs_dst = [];
ker = [];
kIm = [];
unmixing = [];
gmap = [];

for slc=0:SLC-1

    for acc=1:numel(accel)
        R = accel(acc);
        suffix_R = ['n0_s0_slc' num2str(slc) '_encoding_0_R' num2str(R)];
    
        acsSrc = ['GFactor_aug_acsSrc_'  suffix_R];
        acsDst = ['GFactor_aug_acsDst_'  suffix_R];
        ker_str = ['GFactor_aug_ker_'  suffix_R];
        kIm_str = ['GFactor_aug_kIm_'  suffix_R];
        coilMap = ['GFactor_aug_coilMap_'  suffix_R];
        unmixC = ['GFactor_aug_unmixC_'  suffix_R];
        gFactor = ['GFactor_aug_gFactor_'  suffix_R];
    
        acsSrc_array = readGTPlusExportData(fullfile(debug_folder, acsSrc));
        size(acsSrc_array)

        acsDst_array = readGTPlusExportData(fullfile(debug_folder, acsDst));
        size(acsDst_array)

        ker_array = readGTPlusExportData(fullfile(debug_folder, ker_str));
        kIm_array = readGTPlusExportData(fullfile(debug_folder, kIm_str));
        coilMap_array = readGTPlusExportData(fullfile(debug_folder, coilMap));
        unmixC_array = readGTPlusExportData(fullfile(debug_folder, unmixC));
        gFactor_array = niftiread([fullfile(debug_folder, gFactor) '.hdr']);

        acs_src(:,:,:,acc,slc+1) = acsSrc_array;
        acs_dst(:,:,:,acc,slc+1) = acsDst_array;
        %ker(:,:,:,:,acc,slc+1) = ker_array;
        kIm(:,:,:,:,acc,slc+1) = kIm_array;
        unmixing(:,:,:,acc,slc+1) = unmixC_array;
        gmap(:,:,acc,slc+1) = gFactor_array;
    end
end

kspace = single(kspace);
aliasedIm = single(aliasedIm);
acs_src = single(acs_src);
acs_dst = single(acs_dst);
coil_map = single(coil_map);
recon_kspace = single(recon_kspace);
recon_im = single(recon_im);
kIm = single(kIm);
unmixing = single(unmixing);
gmap = single(gmap);