
function [suspected_volume, entropy] = EstimateSuspectedVoxels_Entropy(brainmask, post_csf, post_wm, post_gm, post_outlier, Low_Thres)
% if no one class predominates the segmentation, the voxel is labelled as
% suspected one.

suspected_volume = zeros(size(post_csf), 'uint8');
entropy = zeros(size(post_csf), 'single');

% estimate entropy
index = find( brainmask==1 );
num = length(index);
for i = 1:num
    prob = single([post_csf(index(i)), post_wm(index(i)), post_gm(index(i)), post_outlier(index(i))]+eps) ./ 2048;
    prob = prob + eps;
    if ( sum(find(abs(prob*5-5) <= eps)) >= 1 )
        entropy(index(i)) = 0;
    else
        entropy(index(i)) = -1 * sum(prob.*log2(prob));
    end
end

index = find( ( brainmask==1 ) ...
    & (entropy>=Low_Thres)  );

suspected_volume(index(:)) = 1;

return;