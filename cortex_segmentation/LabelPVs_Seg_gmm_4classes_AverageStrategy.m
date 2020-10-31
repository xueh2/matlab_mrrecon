

function [LabeledSeg, mix] = LabelPVs_Seg_gmm_4classes_AverageStrategy(mix, segResult, ...
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
wmlabel = p;

LabeledSeg = segResult;

ndata = size(mix.priors, 1);
for pp = 1:ndata
    
    j = mix.indexes(pp, 1);
    i = mix.indexes(pp, 2);
    k = mix.indexes(pp, 3);
    
    label = segResult(j, i, k);
    % if the pixel is nonbrain or csf or background (labelled by
    % 0), continue
    if ( (label == 0) | (label == csflabel) | (label == nonbrainlabel) )
        continue;
    end

    [neighbourhood, neighbors] = GetNeighbourhood2(segResult, mix.header, i, j, k, neighborwidth);
%     neighbors = GetNeighbourhood2_6neighbors(segResult, mix.header, mix.offset, i, j, k);
    if ( label == wmlabel )

        % between csf and non-brain tissue
        if ( (isempty(find(neighbors==csflabel)) == 0) & ...
             (isempty(find(neighbors==nonbrainlabel)) == 0) )
            LabeledSeg(j, i, k) = pvlabel;
            % decrease wm/gm priors
            priors = mix.priors(pp, :);
            priors(2) = priors(2)*lamda;
            priors(3) = priors(3)*lamda;
            %priors = normalizePriors(priors);
            temp = 1 - sum(priors);
            %csfNew = priors(1) + temp * priors(1) / (priors(1)+priors(5)+eps);
            csfNew = priors(1) + temp;
            %nonbrainNew = priors(5) + temp * priors(5) / (priors(1)+priors(5)+eps);
            priors(1) = csfNew;
%             priors(5) = nonbrainNew;
            mix.priors(pp, :) = priors;
            
            continue;
        end

        % between csf and gm
        if ( (isempty(find(neighbors==csflabel)) == 0) & ...
             (isempty(find(neighbors==cortexlabel)) == 0) )
            LabeledSeg(j, i, k) = pvlabel;
            % decrease wm priors
            priors = mix.priors(pp, :);
            priors(3) = priors(3)*lamda;
%             priors = normalizePriors(priors);

            % the prior of non-brain tissue can not be incearsed ...
            temp = 1 - sum(priors);
%             csfNew = priors(1) + temp * priors(1) / (priors(1)+priors(2)+eps);
%             gmNew = priors(2) + temp * priors(2) / (priors(1)+priors(2)+eps);
            
            csfNew = priors(1) + temp/2;
            gmNew = priors(2) + temp/2;

            priors(1) = csfNew;
            priors(2) = gmNew;
            mix.priors(pp, :) = priors;
            
            continue;
        end

        % between non-brain tissue and gm
        if ( (isempty(find(neighbors==nonbrainlabel)) == 0) & ...
             (isempty(find(neighbors==cortexlabel)) == 0) )
            LabeledSeg(j, i, k) = pvlabel;
            % decrease wm priors
            priors = mix.priors(pp, :);
            priors(3) = priors(3)*lamda;
%             priors = normalizePriors(priors);
            temp = 1 - sum(priors);
            
%             csfNew = priors(1) + temp * priors(1) / (priors(1)+priors(2)+eps);
%             gmNew = priors(2) + temp * priors(2) / (priors(1)+priors(2)+eps);
            
            csfNew = priors(1) + temp/2;
            gmNew = priors(2) + temp/2;

            priors(1) = csfNew;
            priors(2) = gmNew;
            mix.priors(pp, :) = priors;
            
            continue;
        end

    end

    if ( label == cortexlabel )

        % between csf and non-brain tissue
        if ( (isempty(find(neighbors==csflabel)) == 0) & ...
             (isempty(find(neighbors==nonbrainlabel)) == 0) )
            LabeledSeg(j, i, k) = pvlabel;
            % decrease wm/gm priors
            priors = mix.priors(pp, :);
            priors(2) = priors(2)*lamda;
            priors(3) = priors(3)*lamda;
%             priors = normalizePriors(priors);
            temp = 1 - sum(priors);
%             csfNew = priors(1) + temp * priors(1) / (priors(1)+priors(5)+eps);
            csfNew = priors(1) + temp;
%             nonbrainNew = priors(5) + temp * priors(5) / (priors(1)+priors(5)+eps);
            priors(1) = csfNew;
%             priors(5) = nonbrainNew;
            mix.priors(pp, :) = priors;
            
            continue;
        end

        if ( (isempty(find(neighbors==csflabel)) == 0) & ...
             (isempty(find(neighbors==0)) == 0) )
            LabeledSeg(j, i, k) = pvlabel;
            % decrease gm priors
            priors = mix.priors(pp, :);
            priors(2) = priors(2)*lamda;
%             priors = normalizePriors(priors);
            temp = 1 - sum(priors);
%             csfNew = priors(1) + temp * priors(1) / (priors(1)+priors(5)+eps);
            csfNew = priors(1) + temp;
%             nonbrainNew = priors(5) + temp * priors(5) / (priors(1)+priors(5)+eps);
            priors(1) = csfNew;
            %priors(5) = nonbrainNew;
            mix.priors(pp, :) = priors;
            
            continue;
        end

        if ( (isempty(find(neighbors==wmlabel)) == 0) & ...
             (isempty(find(neighbors==0)) == 0) )
            LabeledSeg(j, i, k) = pvlabel;
            % decrease gm priors
            priors = mix.priors(pp, :);
            priors(2) = priors(2)*lamda;
%             priors = normalizePriors(priors);
%             nonbrainNew = priors(5) + temp * priors(5) / (priors(1)+priors(5));
            temp = 1 - sum(priors);
            priors(1) = priors(1) + temp;
%             priors(5) = nonbrainNew;
            mix.priors(pp, :) = priors;
            
            continue;
        end
    end
end
return
