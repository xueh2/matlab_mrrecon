function [meanRR, meanTriggerTimeDelta, startingPhase, startingInd, endingPhase, endingInd, indexes] = ComputeMeanRR(triggerTime)
% function [meanRR, meanTriggerTimeDelta, startingPhase, endingPhase, indexes] = ComputeMeanRR(triggerTime)
%
% this function detects the cardiac phases in the triggerTime and computes
% the mean R-R interval, mean delta trigger delay
% for every cardiac phase, its staring trigger time is recorded in
% startingPhase and ending trigger time is recorded in endingPhase
% indexes records the image index for every cardiac cycle
%
% Inputs:
%    triggerTime : a group of trigger time continously recorded, including multiple cardiac phases
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

N = numel(triggerTime);

% find the inversion point of trigger time
inversion = [];
for kk=1:N-1    
    if ( triggerTime(kk+1) - triggerTime(kk) > 0 )
        continue;        
    end
    inversion = [inversion kk];    
end
M = numel(inversion);

if ( M < 2 )
    meanRR = 0;
    meanTriggerTimeDelta = 0;
    startingPhase = []; 
    endingPhase = [];
    indexes = [];
end

% fill the starting and ending
startingPhase = zeros(M-1, 1);
endingPhase = zeros(M-1, 1);
startingInd = zeros(M-1, 1);
endingInd = zeros(M-1, 1);

for kk=1:M
    if ( kk ~= M )
        startingPhase(kk) = triggerTime(inversion(kk) + 1);
        startingInd(kk) = inversion(kk) + 1;
    end
    
    if ( kk ~= 1 )
        endingPhase(kk-1) = triggerTime(inversion(kk));
        endingInd(kk-1) = inversion(kk);
    end
end

meanRR = 0;
meanTriggerTimeDelta = [];
indexes = cell(M-1, 1);
for kk=1:M-1
    indexes{kk} = startingInd(kk):endingInd(kk);
    meanRR = meanRR + endingPhase(kk)-startingPhase(kk);
    
    currTriggerTime = triggerTime(indexes{kk});
    meanTriggerTimeDelta = [meanTriggerTimeDelta; currTriggerTime(2:end)-currTriggerTime(1:end-1)];
end
meanRR = meanRR / (M-1);
meanTriggerTimeDelta = mean(meanTriggerTimeDelta);
