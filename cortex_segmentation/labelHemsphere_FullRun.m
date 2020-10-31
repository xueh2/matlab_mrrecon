
function [labelHemsphere, headerLabel, minSlice, maxSlice] = labelHemsphere_FullRun(trans_seg, header, theta, transi)
% the input data is the segmented binary volume in the transverse
% orientation

[coronal, headerCor]=xy2zx(trans_seg, header);
temp = headerCor.xsize;
headerCor.xsize = headerCor.ysize;
headerCor.ysize = temp;
coronal2 = zeros([headerCor.ysize headerCor.xsize headerCor.zsize], 'uint32');
for i = 1:headerCor.zsize
    coronal2(:,:,i) = coronal(:,:,i)';
    
    reflected_slice = reflectSlice_2D(coronal2(:,:,i), 1);
    
    coronal2(:,:,i) = reflected_slice;
end

% for i = 1:headerCor.zsize
%     temp = coronal2(headerCor.ysize+1-i,:,:);
%     coronal2(headerCor.ysize+1-i,:,:) = coronal2(i,:,:);
%     coronal2(i,:,:) = temp;
% end


SaveAnalyze(uint32(coronal2), headerCor, 'wm_seg_roi_Cor.hdr', 'Grey');

% % get the left and right hemisphere (1 and 2)
filename = 'wm_seg_roi_Cor.hdr';
[data, header] = LoadAnalyze(filename, 'Grey');

[labelHemsphere, minSlice, maxSlice] = LabelHemisphere(data, header, theta, transi);
headerLabel = header;
SaveAnalyze(uint32(labelHemsphere), header, 'labelHemsphere.hdr', 'Grey');
return;