function data = load_mrd_output_array(fname)
% data = load_mrd_output_array(fname)

r = mrd.binary.MrdReader(fname)
header = r.read_header();
data = [];
i = 1;
while r.has_data()
    item = r.read_data();
    im = squeeze(item.value);

    if ndims(im)==2
        data(:,:,i) = im;
    end

    if ndims(im)==3
        data(:,:,:,i) = im;
    end

    if ndims(im)==4
        data(:,:,:,:,i) = im;
    end

    i = i + 1;
end
r.close();

disp(['load ' fname ' - ' num2str(size(data))]);