
function LabeledSeg = LabelPVs_Seg_slow_fillLabeledSeg(segResult, header, ...
                                    csflabel, wmlabel, cortexlabel, pvlabel,...
                                    nonbrainlabel, neighborwidth)
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

p = wmlabel(1);
wmlabel = p

LabeledSeg = segResult;

xsize = header.xsize;
ysize = header.ysize;
zsize = header.zsize;

% ndata = size(mix.priors, 1);
for k = 1:zsize
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
                    LabeledSeg(j, i, k) = csflabel;
                    continue;
                end

                % between csf and gm
                if ( (isempty(find(neighbors==csflabel)) == 0) & ...
                     (isempty(find(neighbors==cortexlabel)) == 0) )
                    LabeledSeg(j, i, k) = cortexlabel;
                    continue;
                end

                % between non-brain tissue and gm
                if ( (isempty(find(neighbors==nonbrainlabel)) == 0) & ...
                     (isempty(find(neighbors==cortexlabel)) == 0) )
                    LabeledSeg(j, i, k) = cortexlabel;
                    continue;
                end

            end

            if ( label == cortexlabel )

                % between csf and non-brain tissue
                if ( (isempty(find(neighbors==csflabel)) == 0) & ...
                     (isempty(find(neighbors==nonbrainlabel)) == 0) )
                    LabeledSeg(j, i, k) = csflabel;
                    continue;
                end

                if ( (isempty(find(neighbors==csflabel)) == 0) & ...
                     (isempty(find(neighbors==0)) == 0) )
                    LabeledSeg(j, i, k) = csflabel;
                    continue;
                end

                if ( (isempty(find(neighbors==wmlabel)) == 0) & ...
                     (isempty(find(neighbors==0)) == 0) )
                    LabeledSeg(j, i, k) = wmlabel;
                    continue;
                end
            end
        end
    end
end
