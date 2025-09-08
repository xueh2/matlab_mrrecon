function [data, gmap, res] = prepare_for_stcnnt_inference_spine(debug_dir)
% [data, gmap] = prepare_for_stcnnt_inference_spine(debug_dir)

data = [];
gmap = [];
res = [];

AVE = 15;

for ave=1:AVE
    gmap_name = fullfile(debug_dir, ['input_gmap_stcnnt_set_0_ave_' num2str(ave-1) '_row_0.hdr']);
    if exist(gmap_name)
        gmap(:,:,:,ave) = analyze75read(gmap_name);
        data(:,:,:,ave) = readGTPlusExportData(fullfile(debug_dir, ['input_complex_images_stcnnt_set_0_ave_' num2str(ave-1) '_row_0']) );
        res(:,:,:,ave) = readGTPlusExportData(fullfile(debug_dir, ['output_complex_images_stcnnt_set_0_ave_' num2str(ave-1) '_row_0']) );
    end
end

size(data)
size(gmap)

writeNPY(single(gmap), fullfile(debug_dir,'gmap'));
writeNPY(real(data), fullfile(debug_dir,'input_real'));
writeNPY(imag(data), fullfile(debug_dir,'input_imag'));