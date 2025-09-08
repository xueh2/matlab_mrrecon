function [data, gmap, res] = prepare_for_stcnnt_inference_3D(debug_dir)
% [data, gmap] = prepare_for_stcnnt_inference_3D(debug_dir)

data = [];
gmap = [];
res = [];

SLC = 15;

for slc=1:SLC
    for set=1:2
        gmap_name = fullfile(debug_dir, ['input_gmap_stcnnt_set_' num2str(set-1) '_ave_0_row_' num2str(slc-1) '.hdr']);
        if exist(gmap_name)
            gmap(:,:,:,set, slc) = analyze75read(gmap_name);
            data(:,:,:,set, slc) = readGTPlusExportData(fullfile(debug_dir, ['input_complex_images_stcnnt_set_' num2str(set-1) '_ave_0_row_' num2str(slc-1)]) );
            res(:,:,:,set, slc) = readGTPlusExportData(fullfile(debug_dir, ['output_complex_images_stcnnt_set_' num2str(set-1) '_ave_0_row_' num2str(slc-1)]) );
        end
    end
end

size(data)
size(gmap)

writeNPY(single(gmap), fullfile(debug_dir,'gmap'));
writeNPY(real(data), fullfile(debug_dir,'input_real'));
writeNPY(imag(data), fullfile(debug_dir,'input_imag'));