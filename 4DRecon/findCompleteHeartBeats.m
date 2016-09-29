function [indHB, ind] = findCompleteHeartBeats(ind, indexes)
% find the complete heart beats occupied by some given frames
% ind : the indexes of frames
% indexes: the indexes for every heart beats

numOfHB = numel(indexes);

minHB = 1;
maxHB = numOfHB;

for k=1:numOfHB
    indCurrHB = indexes{k};
    
    if ( ~isempty( find(indCurrHB == ind(1)) ) ) % find the minimal HB
        minHB = k
    end
    
    if ( ~isempty( find(indCurrHB == ind(end)) ) ) % find the minimal HB
        maxHB = k
    end
end

indHB = minHB:maxHB;

ind = [];
for k=minHB:maxHB    
    ind = [ind indexes{k}];
end

