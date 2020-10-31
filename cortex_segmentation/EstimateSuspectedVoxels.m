
function [suspected_volume, maxPosterior] = EstimateSuspectedVoxels(brainmask, post_csf, post_wm, post_gm, post_outlier, Low_Thres)
% if no one class predominates the segmentation, the voxel is labelled as
% suspected one.

suspected_volume = zeros(size(post_csf), 'uint8');
maxPosterior = zeros(size(post_csf), 'uint32');

index = find( ( brainmask==1 ) ...
    & ((post_csf<=Low_Thres) & (post_wm<=Low_Thres) & (post_gm<=Low_Thres) & (post_outlier<=Low_Thres)) );

suspected_volume(index(:)) = 1;

num = length(index);
for i = 1:num
    maxPosterior(index(i)) = max([post_csf(index(i)), post_wm(index(i)), post_gm(index(i)), post_outlier(index(i))]);
end
return;