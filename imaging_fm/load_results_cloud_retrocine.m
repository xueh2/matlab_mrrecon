function [kspace, ref_kspace, kIm, unmixing, coil_map, gmap, im, im_cha, res, header] = load_results_cloud_retrocine(res_dir, patient_id, study_id, meas_id, plot_flag)
% [kspace, ref_kspace, kIm, unmixing, coil_map, gmap, im, im_cha, res, header] = load_results_cloud_retrocine(res_dir, patient_id, study_id, meas_id, plot_flag)
% im, im_cha are complex image, complex multichannel image, res is recon
% final output

kspace = [];
ref_kspace = [];
kIm = [];
unmixing = [];
coil_map = [];
gmap = [];
im = [];
im_cha = [];
res = [];
header = [];

[name, num] = findFILE(res_dir, ['*' patient_id '_' study_id '_' meas_id '*-000000.mrd']);
[res, header] = read_mrd(name{1}, 2048);
if ~isempty(header)
    disp(['load final recon output ' name{1} ' - ' num2str(size(res))]);
    SLC = header.encoding.encoding_limits.slice.maximum + 1;
    PHS = size(res, 3) / SLC;
    res = reshape(res, [size(res,1), size(res,2), PHS, SLC]);
else
    return
end

[name, num] = findFILE(res_dir, ['*' patient_id '_' study_id '_' meas_id '_*coil_map.mrd']);
coil_map = read_mrd(name{1}, 2048);
disp(['load coil map ' name{1} ' - ' num2str(size(coil_map))]);

[name, num] = findFILE(res_dir, ['*' patient_id '_' study_id '_' meas_id '_*gmap.mrd']);
gmap = read_mrd(name{1}, 32);
disp(['load gmap ' name{1} ' - ' num2str(size(gmap))]);

[name, num] = findFILE(res_dir, ['*' patient_id '_' study_id '_' meas_id '_*grappa_unmixing.mrd']);
unmixing = read_mrd(name{1}, 2048);
disp(['load unmixing ' name{1} ' - ' num2str(size(unmixing))]);

[name, num] = findFILE(res_dir, ['*' patient_id '_' study_id '_' meas_id '_*image_domain_kernels.mrd']);
kIm = read_mrd(name{1}, 2048);
disp(['load image_domain_kernels ' name{1} ' - ' num2str(size(kIm))]);

[name, num] = findFILE(res_dir, ['*' patient_id '_' study_id '_' meas_id '_*kspace.mrd']);
kspace = read_mrd(name{1}, 2048);
disp(['load kspace ' name{1} ' - ' num2str(size(kspace))]);

[name, num] = findFILE(res_dir, ['*' patient_id '_' study_id '_' meas_id '_*ref_kspace.mrd']);
ref_kspace = read_mrd(name{1}, 2048);
disp(['load ref_kspace ' name{1} ' - ' num2str(size(ref_kspace))]);

[name, num] = findFILE(res_dir, ['*' patient_id '_' study_id '_' meas_id '_*recon_complex.mrd']);
im = read_mrd(name{1}, 2048);
disp(['load recon_complex ' name{1} ' - ' num2str(size(im))]);

[name, num] = findFILE(res_dir, ['*' patient_id '_' study_id '_' meas_id '_*recon_channel_complex.mrd']);
im_cha = read_mrd(name{1}, 2048);
disp(['load recon_channel_complex ' name{1} ' - ' num2str(size(im_cha))]);

if plot_flag
    s_ratio = 8;
    figure('Name', ['coil map - ' num2str(size(coil_map))]); imagescn(abs(squeeze(coil_map)), [], [], [s_ratio], 3);
    figure('Name', ['gmap - ' num2str(size(gmap))]); imagescn(abs(squeeze(gmap)), [0 2], [], [s_ratio]);
    figure('Name', ['unmixing - ' num2str(size(unmixing))]); imagescn(abs(squeeze(unmixing)), [], [], [s_ratio], 4);
    figure('Name', ['kIm - ' num2str(size(kIm))]); imagescn(abs(squeeze(kIm(:,:,:,1,1,1,:))), [], [], [s_ratio], 4);
    figure('Name', ['kspace - ' num2str(size(kspace))]); imagescn(abs(squeeze(kspace(:,:,1,1,1,1,:))), [], [], s_ratio);PerfColorMap;
    figure('Name', ['ref_kspace - ' num2str(size(ref_kspace))]); imagescn(abs(squeeze(ref_kspace(:,:,1,1,1,1,:))), [], [], s_ratio);colormap('jet');
    figure('Name', ['im - ' num2str(size(im))]); imagescn(abs(squeeze(im)), [], [], [s_ratio], 3);
    figure('Name', ['im_cha - ' num2str(size(im_cha))]); imagescn(abs(squeeze(im_cha(:,:,:,1,1,:))), [], [], [s_ratio], 3);
    if ~isempty(res)
        figure('Name', ['res - ' num2str(size(res))]); imagescn(abs(squeeze(res)), [], [], [s_ratio], 3);
    end
end

coil_map = squeeze(coil_map);
gmap = squeeze(gmap);
unmixing = squeeze(unmixing);
kIm = squeeze(kIm);
kspace = squeeze(kspace);
ref_kspace = squeeze(ref_kspace);
im = squeeze(im);
im_cha = squeeze(im_cha);
if ~isempty(res)
    res = squeeze(res);
end

end

function [data, header] = read_mrd(fname, min_size)

    data = [];
    header = [];

    s = dir(fname);         
    filesize = s.bytes;

    if filesize < min_size*1024
        return
    end

    r = mrd.binary.MrdReader(fname);
    header = r.read_header();
    data = [];
    i = 1;
    while r.has_data()
        item = r.read_data();
        try
            im = item.value.data;
        catch
            im = item.value;
        end
        if ndims(im) == 2
            data(:,:,i) = im;
        end
        if ndims(im) == 3
            data(:,:,:,i) = im;
        end
        if ndims(im) == 4
            data(:,:,:,:,i) = im;
        end
        if ndims(im) == 5
            data(:,:,:,:,:,i) = im;
        end

        if ndims(im) == 6
            data(:,:,:,:,:,:,i) = im;
        end

        if ndims(im) == 7
            data(:,:,:,:,:,:,:,i) = im;
        end
        if ndims(im) == 8
            data(:,:,:,:,:,:,:,:,i) = im;
        end

        i = i + 1;
    end
    r.close();
end