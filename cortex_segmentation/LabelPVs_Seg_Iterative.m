
function [LabeledSeg, segResult] = LabelPVs_Seg_Iterative(segResult, header, ...
                                    csflabel, wmlabel, cortexlabel, pvlabel,...
                                    nonbrainlabel, neighborwidth, maxIter)
% according to the prior knowledge to label the pvs from the initial
% segmentation
% the pvs arel labelled as pvlabel

disp('Label partial volumes ...');

segResult2 = segResult;
num = length(wmlabel);
if (num>1)
    
    for i = 2:num
        segResult2(find(segResult==wmlabel(i))) = wmlabel(1);
    end
    
end
wmlabelOld = wmlabel;
p = wmlabel(1);
wmlabel = p

LabeledSeg = segResult;

xsize = header.xsize;
ysize = header.ysize;
zsize = header.zsize;

changed = true;
numIter = 1;
while ( changed )
    disp(['The iteration ' num2str(numIter) ' ... ']);
    if ( numIter > maxIter )
        break;
    end
    changed = false;
    for k = 1:zsize
        k
        for j = 1:ysize
            for i = 1:xsize

                label = segResult2(j, i, k);
                % if the pixel is nonbrain or csf or background (labelled by
                % 0), continue
                if ( (label == 0) | (label == csflabel) | (label == nonbrainlabel) )
                    continue;
                end

                [neighbourhood, neighbors] = GetNeighbourhood2(segResult2, header, i, j, k, neighborwidth);

                if ( label == wmlabel )

                    % between csf and non-brain tissue
                    if ( (isempty(find(neighbors==csflabel)) == 0) & ...
                         (isempty(find(neighbors==nonbrainlabel)) == 0) )
                        LabeledSeg(j, i, k) = pvlabel;
                        segResult(j, i, k) = csflabel;
                        changed = true;
                    end

                    % between csf and gm
                    if ( (isempty(find(neighbors==csflabel)) == 0) & ...
                         (isempty(find(neighbors==cortexlabel)) == 0) )
                        LabeledSeg(j, i, k) = pvlabel;
                        segResult(j, i, k) = cortexlabel;
                        changed = true;
                    end

                    % between non-brain tissue and gm
                    if ( (isempty(find(neighbors==nonbrainlabel)) == 0) & ...
                         (isempty(find(neighbors==cortexlabel)) == 0) )
                        LabeledSeg(j, i, k) = pvlabel;
                        segResult(j, i, k) = nonbrainlabel;
                        changed = true;
                    end

                end

                if ( label == cortexlabel )

                    % between csf and non-brain tissue
                    if ( (isempty(find(neighbors==csflabel)) == 0) & ...
                         (isempty(find(neighbors==nonbrainlabel)) == 0) )
                        LabeledSeg(j, i, k) = pvlabel;
                        segResult(j, i, k) = csflabel;
                        changed = true;
                    end

                    if ( (isempty(find(neighbors==csflabel)) == 0) & ...
                         (isempty(find(neighbors==0)) == 0) )
                        LabeledSeg(j, i, k) = pvlabel;
                        segResult(j, i, k) = csflabel;
                        changed = true;
                    end

                    if ( (isempty(find(neighbors==wmlabel)) == 0) & ...
                         (isempty(find(neighbors==0)) == 0) )
                        LabeledSeg(j, i, k) = pvlabel;
                        segResult(j, i, k) = 0;
                        changed = true;
                    end
                end
            end
        end
    end
    numIter = numIter + 1;
    
    segResult2 = segResult;
    num = length(wmlabelOld);
    if (num>1)

        for i = 2:num
            segResult2(find(segResult==wmlabelOld(i))) = wmlabelOld(1);
        end

    end
end