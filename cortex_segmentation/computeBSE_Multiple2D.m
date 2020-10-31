
function [BSE, distPercentages, dist] = computeBSE_Multiple2D(seg_result, cortexlabel, ...
            transverseLabel, labeled_trans, ...
            sagittalLabel, labeled_sag, ...
            coronaryLabel, labeled_cor, ... 
            header, distThresholds)
        
% compute the Dice Similarity Measure
% the unused orientation is marked by -1 as its label
% compute the False Positive and False Negative

segLabel = false(size(seg_result));
num = length(cortexlabel);
for i = 1:num
    segLabel(find(seg_result==cortexlabel(i))) = 1;
end

dist = [];
distRealNumbers = zeros(size(distThresholds));

xvoxelsize = header.xvoxelsize;
yvoxelsize = header.yvoxelsize;
zvoxelsize = header.zvoxelsize;

% transverse
if ( isempty(find(labeled_trans==-1)) == 1 )
    num = length(labeled_trans);
    for i = 1:num
        
        auto2DSlice = segLabel(:,:, labeled_trans(i));
        manual2DSlice = transverseLabel(:,:, labeled_trans(i));
        
        [BSE, distPercentages, currentdistRealNumbers, currentdist] = computeBSE_2D(auto2DSlice, manual2DSlice, xvoxelsize, yvoxelsize, distThresholds);
        
        dist = [dist; currentdist];
        distRealNumbers = distRealNumbers + currentdistRealNumbers;
    end
end

% sagittal
if ( isempty(find(labeled_sag==-1)) == 1 )
    num = length(labeled_sag);
    for i = 1:num
        
        auto2DSlice = segLabel(:,labeled_sag(i), :);
        manual2DSlice = sagittalLabel(:,labeled_sag(i), :);

        [BSE, distPercentages, currentdistRealNumbers, currentdist] = computeBSE_2D(auto2DSlice, manual2DSlice, xvoxelsize, zvoxelsize, distThresholds);
        
        dist = [dist; currentdist];
        distRealNumbers = distRealNumbers + currentdistRealNumbers;
    end
end

% coronal
if ( isempty(find(labeled_cor==-1)) == 1 )
    num = length(labeled_cor);
    for i = 1:num
        
        auto2DSlice = segLabel(labeled_cor(i), :, :);
        manual2DSlice = coronaryLabel(labeled_cor(i), :, :);

        [BSE, distPercentages, currentdistRealNumbers, currentdist] = computeBSE_2D(auto2DSlice, manual2DSlice, yvoxelsize, zvoxelsize, distThresholds);
        
        dist = [dist; currentdist];
        distRealNumbers = distRealNumbers + currentdistRealNumbers;
    end
end

BSE = [mean(dist) std(dist) max(dist)];

distThresholds = distThresholds .* (xvoxelsize+yvoxelsize+zvoxelsize)/3;

distPercentages = zeros(size(distThresholds));
distRealNumbers = zeros(size(distThresholds));

absdist = abs(dist);
numdist = length(dist);

for i=1:length(distThresholds)
    distRealNumbers(i) = length( find( absdist>distThresholds(i) ) );
    distPercentages(i) = distRealNumbers(i) / numdist;
end

return;