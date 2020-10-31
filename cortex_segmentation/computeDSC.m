
function dsc = computeDSC(seg_result, cortexlabel, ...
            transverseLabel, labeled_trans, ...
            sagittalLabel, labeled_sag, ...
            coronaryLabel, labeled_cor, ... 
            header)
        
% compute the Dice Similarity Measure
% the unused orientation is marked by -1 as its label


segLabel = false(size(seg_result));
num = length(cortexlabel);
for i = 1:num
    segLabel(find(seg_result==cortexlabel(i))) = 1;
end

segLabelSlices = false(size(seg_result));

if ( isempty(find(labeled_trans==-1)) == 1 )
    num = length(labeled_trans);
    for i = 1:num
        segLabelSlices(:,:, labeled_trans(i)) = segLabel(:,:,labeled_trans(i));
    end
end

if ( isempty(find(labeled_sag==-1)) == 1 )
    num = length(labeled_sag);
    for i = 1:num
        segLabelSlices(:, labeled_sag(i), :) = segLabel(:,labeled_sag(i), :);
    end
end

if ( isempty(find(labeled_cor==-1)) == 1 )
    num = length(labeled_cor);
    for i = 1:num
        segLabelSlices(labeled_cor(i), :, :) = segLabel(labeled_cor(i), :, :);
    end
end
%SaveAnalyze(uint32(segLabelSlices), header, 'labeled.hdr', 'Grey');

segnum = length(find(segLabelSlices==1));

if ( isempty(transverseLabel) == 1 )
    transverseLabel = zeros(size(seg_result), 'uint32');
end
if ( isempty(sagittalLabel) == 1 )
    sagittalLabel = zeros(size(seg_result), 'uint32');
end
if ( isempty(coronaryLabel) == 1 )
    coronaryLabel = zeros(size(seg_result), 'uint32');
end
manualLabel = transverseLabel | sagittalLabel | coronaryLabel;

manualnum = length(find(manualLabel==1));

crossLabel = segLabel & manualLabel;
crossnum = length(find(crossLabel==1));

dsc = 2 * crossnum / (segnum+manualnum);

return;