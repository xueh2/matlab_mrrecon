function [data, gmap] = prepare_for_stcnnt_inference_perf(debug_dir)
% [data, gmap] = prepare_for_stcnnt_inference_perf(debug_dir)

data = [];
gmap = [];

for slc=1:12
    gmap_name = fullfile(debug_dir, ['input_gmap_stcnnt_slc_' num2str(slc-1) '_set_0_row_' num2str(slc-1) '.hdr']);
    if exist(gmap_name)
        gmap(:,:,slc) = analyze75read(gmap_name);
        data(:,:,:,slc) = readGTPlusExportData(fullfile(debug_dir, ['input_complex_images_stcnnt_slc_' num2str(slc-1) '_set_0_row_' num2str(slc-1)]) );
    end
end

size(data)
size(gmap)

writeNPY(single(gmap), fullfile(debug_dir,'gmap'));
writeNPY(real(data), fullfile(debug_dir,'input_real'));
writeNPY(imag(data), fullfile(debug_dir,'input_imag'));