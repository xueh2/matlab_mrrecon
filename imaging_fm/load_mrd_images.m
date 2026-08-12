function [data, all_heads, all_meta, header] = load_mrd_images(fname, series_to_load)
% [data, all_heads, all_meta] = load_mrd_images(fname, series_to_load)

if nargin < 2
    series_to_load = 0;
end

r = mrd.binary.MrdReader(fname)
header = r.read_header();
data = [];
all_heads = [];
all_meta = [];
i = 1;
while r.has_data()
    item = r.read_data();
    im = squeeze(item.value.data);
    h = item.value.head;
    meta = item.value.meta;

    if h.image_series_index ~= series_to_load
        continue;
    end

    if ndims(im)==2
        data(:,:,i) = im;
    end

    if ndims(im)==3
        data(:,:,:,i) = im;
    end

    if ndims(im)==4
        data(:,:,:,:,i) = im;
    end

    all_heads = [all_heads; h];
    all_meta = [all_meta; meta];

    i = i + 1;
end
r.close();

disp(['load ' fname ' - ' num2str(size(data))]);