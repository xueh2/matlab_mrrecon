
% read original and binning kspace
function [ori, binned] = load_binning_debug(debug_dir)

ori = [];
binned = [];

slc = 0;

while(isFileExist(fullfile(debug_dir, ['full_kspace_raw_encoding_0_SLC_0_SLCOrder_' num2str(slc) '_MAG.hdr'])))
    ori{slc+1} = readGTPlusExportData(fullfile(debug_dir, ['full_kspace_raw_encoding_0_SLC_0_SLCOrder_' num2str(slc)]));
    slc = slc + 1;
end

slc = 0;
while(isFileExist(fullfile(debug_dir, ['kspace_binning_S_0_encoding_0_SLC_0_SLCOrder_' num2str(slc) '_MAG.hdr'])))
    binned{slc+1} = readGTPlusExportData(fullfile(debug_dir, ['kspace_binning_S_0_encoding_0_SLC_0_SLCOrder_' num2str(slc)]));
    slc = slc + 1;
end