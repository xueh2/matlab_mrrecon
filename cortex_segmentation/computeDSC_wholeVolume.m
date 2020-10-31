
function dsc = computeDSC_wholeVolume(seg_result, cortexlabel, ...
            manualSeg, mlabel)
        
% compute the Dice Similarity Measure
% the unused orientation is marked by -1 as its label


segLabel = false(size(seg_result));
num = length(cortexlabel);
for i = 1:num
    segLabel(find(seg_result==cortexlabel(i))) = 1;
end

manualLabel = false(size(manualSeg));
num = length(mlabel);
for i = 1:num
    manualLabel(find(manualSeg==mlabel(i))) = 1;
end

segnum = length(find(segLabel==1));
manualnum = length(find(manualLabel==1));

crossLabel = segLabel & manualLabel;
crossnum = length(find(crossLabel==1));

dsc = 2 * crossnum / (segnum+manualnum);

return;