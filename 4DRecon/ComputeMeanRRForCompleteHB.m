function [meanRR, meanTriggerTimeDelta, startingPhase, startingInd, endingPhase, endingInd, indexesHB] = ComputeMeanRRForCompleteHB(indHB, indexes, triggerTime)
% function [meanRR, meanTriggerTimeDelta, startingPhase, startingInd, endingPhase, endingInd, indexes] = ComputeMeanRRForCompleteHB(indHB, indexes)
%
% this function takes the input of complete heart beats and computes
% the mean R-R interval, mean delta trigger delay
% for every cardiac phase, its staring trigger time is recorded in
% startingPhase and ending trigger time is recorded in endingPhase
% indexes records the image index for every cardiac cycle
%
% Inputs:
%    indHB : complete heart beats
%    indexes: the frame indexes for every heart beats
%
% Output:
%    meanRR : mean R-R inverval, computed as the mean duration of all recorded cardiac cycles
%    meanTriggerTimeDelta : trigger time delta is the interval between two consecutive trigger time
%    startingPhase : for every complete cardiac phase included in the input, its starting trigger time is recorded
%    startingInd : index of every starting phase
%    endingPhase : for every complete cardiac phase included in the input, its ending trigger time is recorded
%    endingInd : index of every ending phase
%    indexes : for every complete cardiac phase included in the input, the
%    indexes of all images in this phase is recorded as one item of indexes

%     ***************************************
%     *  Hui Xue (hui-xue@siemens.com       *
%     *  2012-03                            *
%     ***************************************

N = numel(indHB);

% fill the starting and ending
startingPhase = zeros(N, 1);
endingPhase = zeros(N, 1);
startingInd = zeros(N, 1);
endingInd = zeros(N, 1);

firstInd = indexes{indHB(1)};
firstInd = firstInd(1);

for kk=1:N
    ind = indexes{indHB(kk)};    
    startingPhase(kk) = triggerTime(ind(1));
    startingInd(kk) = ind(1);
    endingPhase(kk) = triggerTime(ind(end));
    endingInd(kk) = ind(end);
end

meanRR = 0;
meanTriggerTimeDelta = [];
indexesHB = cell(N, 1);
for kk=1:N
    indexesHB{kk} = startingInd(kk):endingInd(kk);
    meanRR = meanRR + endingPhase(kk)-startingPhase(kk);
    
    currTriggerTime = triggerTime(indexesHB{kk});
    meanTriggerTimeDelta = [meanTriggerTimeDelta; currTriggerTime(2:end)-currTriggerTime(1:end-1)];
end
meanRR = meanRR / N;
meanTriggerTimeDelta = mean(meanTriggerTimeDelta);
