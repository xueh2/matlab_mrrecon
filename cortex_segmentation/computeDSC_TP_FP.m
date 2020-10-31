
function [DSC, TP, FP] = computeDSC_TP_FP(seg_result, manualSeg)
% compute the dice value, Ture Positive and False Positive

sizeSeg = size(seg_result);
sizeManual = size(manualSeg);

if ( same(sizeSeg, sizeManual)==0 )
    error('the size of two inputs should be the same ...');
end

segnum = length(find(seg_result>0));
manualnum = length(find(manualSeg>0));

crossLabel = seg_result & manualSeg;
crossnum = length(find(crossLabel>0));

dsc = 2 * crossnum / (segnum+manualnum);

TP = crossnum / manualnum;
FP = (segnum-crossnum) / manualnum;

return