function [data, gmap, res] = prepare_for_stcnnt_inference_LGE(debug_dir)
% [data, gmap] = prepare_for_stcnnt_inference_LGE(debug_dir)

data = [];
gmap = [];
res = [];

SLC = 23;

% for slc=1:SLC
%     gmap_name = fullfile(debug_dir, ['input_gmap_stcnnt_slc_' num2str(slc-1) '_set_0_con_0_row_' num2str(slc-1) '.hdr']);
%     if exist(gmap_name)
%         gmap(:,:,slc) = analyze75read(gmap_name);
%         data(:,:,:,1, slc) = readGTPlusExportData(fullfile(debug_dir, ['input_complex_images_stcnnt_slc_' num2str(slc-1) '_set_0_con_0_row_' num2str(slc-1)]) );
%         res(:,:,:,1, slc) = readGTPlusExportData(fullfile(debug_dir, ['output_complex_images_stcnnt_slc_' num2str(slc-1) '_set_0_con_0_row_' num2str(slc-1)]) );
%         data(:,:,:,2,slc) = readGTPlusExportData(fullfile(debug_dir, ['input_complex_images_stcnnt_slc_' num2str(slc-1) '_set_1_con_0_row_' num2str(slc-1)]) );
%         res(:,:,:,2,slc) = readGTPlusExportData(fullfile(debug_dir, ['output_complex_images_stcnnt_slc_' num2str(slc-1) '_set_1_con_0_row_' num2str(slc-1)]) );
%     end
% end

for slc=1:SLC
    gmap_name = fullfile(debug_dir, ['input_gmap_stcnnt_slc_' num2str(slc-1) '_set_0_con_0_row_0.hdr']);
    if exist(gmap_name)
        gmap(:,:,slc) = analyze75read(gmap_name);
        data(:,:,:,1, slc) = readGTPlusExportData(fullfile(debug_dir, ['input_complex_images_stcnnt_slc_' num2str(slc-1) '_set_0_con_0_row_0']) );
        res(:,:,:,1, slc) = readGTPlusExportData(fullfile(debug_dir, ['output_complex_images_stcnnt_slc_' num2str(slc-1) '_set_0_con_0_row_0']) );
        data(:,:,:,2,slc) = readGTPlusExportData(fullfile(debug_dir, ['input_complex_images_stcnnt_slc_' num2str(slc-1) '_set_1_con_0_row_0']) );
        res(:,:,:,2,slc) = readGTPlusExportData(fullfile(debug_dir, ['output_complex_images_stcnnt_slc_' num2str(slc-1) '_set_1_con_0_row_0']) );
    end
end

if size(gmap,3) == 1
    gmap = repmat(gmap, [1 1 2]);
end


size(data)
size(gmap)

writeNPY(single(gmap), fullfile(debug_dir,'gmap'));
writeNPY(real(data), fullfile(debug_dir,'input_real'));
writeNPY(imag(data), fullfile(debug_dir,'input_imag'));