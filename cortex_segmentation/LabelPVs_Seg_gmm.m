
function [LabeledSeg, mix] = LabelPVs_Seg_gmm(mix, segResult, header, ...
                                    csflabel, wmlabel, cortexlabel, pvlabel,...
                                    nonbrainlabel)
% according to the prior knowledge to label the pvs from the initial
% segmentation
% the pvs arel labelled as pvlabel

disp('Label partial volumes ...');

neighborwidth = mix.neighborwidth;
lamda = mix.lamda;

num = length(wmlabel);
if (num>1)
    
    for i = 2:num
        segResult(find(segResult==wmlabel(i))) = wmlabel(1);
    end
    
end

p = wmlabel(1);
wmlabel = p

LabeledSeg = segResult;

xsize = header.xsize;
ysize = header.ysize;
zsize = header.zsize;

ndata = size(mix.priors, 1);
for pp = 1:ndata
    
    [j, i, k] = mix.indexes(pp, :);
    
    label = segResult(j, i, k);
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
            % decrease wm/gm priors
            
        end

        % between csf and gm
        if ( (isempty(find(neighbors==csflabel)) == 0) & ...
             (isempty(find(neighbors==cortexlabel)) == 0) )
            LabeledSeg(j, i, k) = pvlabel;
            % decrease wm priors
        end

        % between non-brain tissue and gm
        if ( (isempty(find(neighbors==nonbrainlabel)) == 0) & ...
             (isempty(find(neighbors==cortexlabel)) == 0) )
            LabeledSeg(j, i, k) = pvlabel;
            % decrease wm priors
        end

    end

    if ( label == cortexlabel )

        % between csf and non-brain tissue
        if ( (isempty(find(neighbors==csflabel)) == 0) & ...
             (isempty(find(neighbors==nonbrainlabel)) == 0) )
            LabeledSeg(j, i, k) = pvlabel;
        end

        if ( (isempty(find(neighbors==csflabel)) == 0) & ...
             (isempty(find(neighbors==0)) == 0) )
            LabeledSeg(j, i, k) = pvlabel;
        end

        if ( (isempty(find(neighbors==wmlabel)) == 0) & ...
             (isempty(find(neighbors==0)) == 0) )
            LabeledSeg(j, i, k) = pvlabel;
        end
    end
