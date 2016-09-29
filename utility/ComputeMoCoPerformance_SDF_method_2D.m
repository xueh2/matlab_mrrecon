function [frame_index, dice, FP, FN, BSE] = ComputeMoCoPerformance_SDF_method_2D(seg2D, segSDF_moco, voxelSize)
% Compute the moco performance statistics
% frame_index: record which frames have the segmented contours
% dice: M*M dice ratio matrix, same as FP, FN and BSE
% BSE: M*M*3, for mean, std, and max

N = size(seg2D, 3);

frame_index = [];
for n=1:N
    s2D = seg2D(:,:,n);
    s2D = double(s2D);
    
    ind = find(s2D(:)>0);
    if (~isempty(ind))
        frame_index = [frame_index; n];
    end
end

num = numel(frame_index);

dice = ones(num, num);
FP = zeros(num, num);
FN = zeros(num, num);
BSE = zeros(num, num, 3);

for n=1:num
    n
    for m=n+1:num
        
        sdf1 = segSDF_moco(:,:,frame_index(n));
        sdf2 = segSDF_moco(:,:,frame_index(m));
        
        s1 = zeros(size(sdf1));
        s2 = zeros(size(sdf2));
        
        s1(find(sdf1<0)) = 1;
        s2(find(sdf2<0)) = 1;
        
        [a_dice, a_FP, a_FN] = ComputeDSC_2D_FP_FN(s1, s2);
        [a_BSE, distPercentages, distRealNumbers, dist] = computeBSE_WithSDF_2D(s1, sdf1, s2, sdf2, voxelSize(1), voxelSize(2), 1.5*voxelSize(1));
        a_BSE = 0;
        
        dice(n, m) = a_dice;
        FP(n, m) = a_FP;
        FN(n, m) = a_FN;
        BSE(n, m, :) = a_BSE(:);
        
        dice(m, n) = a_dice;
        FP(m, n) = a_FP;
        FN(m, n) = a_FN;
        BSE(m, n, :) = a_BSE(:);
    end
end
