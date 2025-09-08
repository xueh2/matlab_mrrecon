function [data, gmap, res] = prepare_for_stcnnt_inference_FW(debug_dir)
% [data, gmap] = prepare_for_stcnnt_inference_LGE(debug_dir)

data = [];
gmap = [];
res = [];

SLC = 3;
CON = 4;

for slc=1:SLC
    for con=1:CON
        gmap_name = fullfile(debug_dir, ['input_gmap_stcnnt_slc_' num2str(slc-1) '_set_0_con_' num2str(con-1) '_row_0.hdr']);
        if exist(gmap_name)        
            if exist(fullfile(debug_dir, ['output_complex_images_stcnnt_slc_' num2str(slc-1) '_set_0_con_' num2str(con-1) '_row_0_MAG.hdr']))
                gmap(:,:,con, slc) = analyze75read(gmap_name);
                dd = readGTPlusExportData(fullfile(debug_dir, ['input_complex_images_stcnnt_slc_' num2str(slc-1) '_set_0_con_' num2str(con-1) '_row_0']) );;
                if numel(size(dd))==3
                    data(:,:,:,slc) = dd;
                    res(:,:,:, slc) = readGTPlusExportData(fullfile(debug_dir, ['output_complex_images_stcnnt_slc_' num2str(slc-1) '_set_0_con_' num2str(con-1) '_row_0']) );
                else
                    data(:,:,:,con, slc) = readGTPlusExportData(fullfile(debug_dir, ['input_complex_images_stcnnt_slc_' num2str(slc-1) '_set_0_con_' num2str(con-1) '_row_0']) );
                    res(:,:,:,con, slc) = readGTPlusExportData(fullfile(debug_dir, ['output_complex_images_stcnnt_slc_' num2str(slc-1) '_set_0_con_' num2str(con-1) '_row_0']) );
                end
            end
        end
    end
end

size(data)
size(gmap)

writeNPY(single(gmap), fullfile(debug_dir,'gmap'));
writeNPY(real(data), fullfile(debug_dir,'input_real'));
writeNPY(imag(data), fullfile(debug_dir,'input_imag'));