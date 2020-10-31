
function newGlobal_SegResult = Combine_Lodal_Global(global_SegResult, tissuelabels_seg,...
    partsNumber, local_SegResults)
% combine the local segmentation into the global results

numTissue = length(tissuelabels_seg);
newGlobal_SegResult = global_SegResult;

for k = 1:numTissue
    newGlobal_SegResult(find(newGlobal_SegResult == tissuelabels_seg(k))) = 0;
end

for i = 1:partsNumber
    i
    for k = 1:numTissue
        newGlobal_SegResult(find( local_SegResults{i} == tissuelabels_seg(k))) = tissuelabels_seg(k);
    end
end

return;