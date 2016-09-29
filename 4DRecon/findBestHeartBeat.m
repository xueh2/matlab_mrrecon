function [bestCardiacPhase, costFunc] = findBestHeartBeat(cardiacPhases, meanHBSlcLocation, binData, indTrigger, triggerTimeSLC, ... 
    prevHBImAtRequiredCardiacPhases, prevBestHBSlcLocation, slcLocationIncrease, strategy, mask)
% find the best heart beat, given the heart beat images from previous cardiac cycle
% Four different strategies are implemented
%
% to avoid the influences of floating slice location, a mask can be supplied to cover only the heart and diaphragm region
%
% 1) 'SSD': SSD + slice shifting constraint: the heart beat with minimal SSD to the previous heart beat is picked
%
% 2) 'MOCO': moco + slice shifting constraint: the heart beat with minimal total deformation to the prev heart beat is picked
%
% 3) 'MOCO_NoSliceShifting': moco without slice shifting constraint: the heart beat with minimal total deformation to the prev heart beat is picked; the mean slice location
% of the selected heart beat is not required to be upper or lower than previous best heart beat
%
% The slice shifting constraint: for every current heart beat, first compare its mean slice location to the prev best heart beat 
% if the slice location moves in the correct direction, this heart beat is lised as candidates

numOfCardiacCycleInSLC = numel(indTrigger);
Nfe = size(prevHBImAtRequiredCardiacPhases, 1);
Npe = size(prevHBImAtRequiredCardiacPhases, 2);

mask = repmat(mask, [1 1 numel(cardiacPhases)]);
indUsed = find(mask(:)>0);

if ( strcmp(lower(strategy), 'ssd')==1 )
    
    % compute the mean of total SSD between every current heart beat and previous best heart beat
    meanTotalSSD = zeros(numOfCardiacCycleInSLC, 1);
    meanTotalSSD(:) = -1;

    % First, compute the interpolated images of the current heart beat at the required cardiac phases
    % then, compute the mean SSD between each pair of corresponding cardiac phases
    % third, apply the mask and compute the SSD

    % for every current HB
    firstValideHB = 0;
    for b=1:numOfCardiacCycleInSLC

        if ( slcLocationIncrease )
            if ( prevBestHBSlcLocation >= meanHBSlcLocation(b) )
                continue;
            end
        else
            if ( prevBestHBSlcLocation <= meanHBSlcLocation(b) )
                continue;
            end
        end

        if ( firstValideHB==0 )
            firstValideHB = b;
        end
        
        indForCurrHB = indTrigger{b};
        binDataForCurrHB = binData(:,:,indForCurrHB);
        triggerTimeForCurrHB = triggerTimeSLC(indForCurrHB);

        currHBImAtRequiredCardiacPhases = performCineInterpolatiorTriggerTime(binDataForCurrHB, triggerTimeForCurrHB, cardiacPhases);

        diffHB = prevHBImAtRequiredCardiacPhases(indUsed(:))-currHBImAtRequiredCardiacPhases(indUsed(:));
        % diffHB = prevHBImAtRequiredCardiacPhases-currHBImAtRequiredCardiacPhases;
        meanTotalSSD(b) = norm(diffHB(:));
    end

    bestCardiacPhase = -1;
    minTotalSSD = inf;
    for b=1:numOfCardiacCycleInSLC
        if ( meanTotalSSD(b) ~= -1 )
            if ( meanTotalSSD(b) < minTotalSSD )
                bestCardiacPhase = b;
                minTotalSSD = meanTotalSSD(b);
            end
        end
    end
    
    % the next heart beat should not be two far way from the current one
    if ( (bestCardiacPhase-firstValideHB) >= 10 )
        [meanTotalSSD, bestCardiacPhase] = min(meanTotalSSD(firstValideHB:firstValideHB+10));
        bestCardiacPhase = bestCardiacPhase + firstValideHB - 1;
    end
    
    costFunc = meanTotalSSD;
    
elseif ( strcmp(lower(strategy), 'moco')==1 )
    
    % moco parameters
    inverse = 1;
    initial = 0;
    iters = [64 64 8];
    sigma = Nfe/16+Npe/16;
    neighbor = 2.0;
    stepDiv = 3.0;
    moreIterInv = 1;
    algo = 'GLCC';
    volumePreserving = 0;
        
    meanTotalDeform = zeros(numOfCardiacCycleInSLC, 1);
    meanTotalDeform(:) = inf;
    
    % First, compute the interpolated images of the current heart beat at the required cardiac phases
    % then, perform the moco between corresponding cardiac phase images between current heart beat and previous heart beat
    % third, apply the mask and compute the mean deformation
    
    % for every current HB
    firstValideHB = 0;
    for b=1:numOfCardiacCycleInSLC

        if ( slcLocationIncrease )
            if ( prevBestHBSlcLocation >= meanHBSlcLocation(b) )
                continue;
            end
        else
            if ( prevBestHBSlcLocation <= meanHBSlcLocation(b) )
                continue;
            end
        end

        if ( firstValideHB==0 )
            firstValideHB = b;
        end
        
        indForCurrHB = indTrigger{b};
        binDataForCurrHB = binData(:,:,indForCurrHB);
        triggerTimeForCurrHB = triggerTimeSLC(indForCurrHB);

        currHBImAtRequiredCardiacPhases = performCineInterpolatiorTriggerTime(binDataForCurrHB, triggerTimeForCurrHB, cardiacPhases);

        % perform the moco
        header = CreateFtkHeaderInfo(currHBImAtRequiredCardiacPhases, [1.8 1.8 3.6]);
        deformInitial = zeros(size(currHBImAtRequiredCardiacPhases));
        
        [moco, dx, dy, dxInv, dyInc] = Matlab_PerformTemporalPairWiseMotionCorrection(double(prevHBImAtRequiredCardiacPhases), double(currHBImAtRequiredCardiacPhases), ...
            header, inverse, initial, iters, sigma, neighbor, stepDiv, moreIterInv, algo, deformInitial, deformInitial, deformInitial, deformInitial, volumePreserving);
        
        dxSqr = dx.*dx;
        dySqr = dy.*dy;
        sumDeform = dxSqr(indUsed(:)) + dySqr(indUsed(:));
        % sumDeform = dxSqr(:) + dySqr(:);
        meanTotalDeform(b) = norm(sumDeform);
    end

    bestCardiacPhase = -1;
    minTotalDeform = inf;
    for b=1:numOfCardiacCycleInSLC
        if ( meanTotalDeform(b) ~= inf )
            if ( meanTotalDeform(b) < minTotalDeform )
                bestCardiacPhase = b;
                minTotalDeform = meanTotalDeform(b);
            end
        end
    end    
    
    % [minTotalDeform, bestCardiacPhase] = min(meanTotalDeform);
    
    % the next heart beat should not be two far way from the current one
    if ( (bestCardiacPhase-firstValideHB) >= 10 )
        [minTotalDeform, bestCardiacPhase] = min(meanTotalDeform(firstValideHB:firstValideHB+10));
        bestCardiacPhase = bestCardiacPhase + firstValideHB - 1;
    end
    
    costFunc = meanTotalDeform;
    
elseif ( strcmp(lower(strategy), 'moco_nosliceshifting')==1 )
    
    % moco parameters
    inverse = 1;
    initial = 0;
    iters = [64 64 0];
    sigma = Nfe/16+Npe/16;
    neighbor = 2.0;
    stepDiv = 3.0;
    moreIterInv = 1;
    algo = 'GLCC';
    volumePreserving = 0;
        
    meanTotalDeform = zeros(numOfCardiacCycleInSLC, 1);
    meanTotalDeform(:) = inf;
    
    % First, compute the interpolated images of the current heart beat at the required cardiac phases
    % then, perform the moco between corresponding cardiac phase images between current heart beat and previous heart beat
    % third, apply the mask and compute the mean deformation
    
    % for every current HB
    for b=1:numOfCardiacCycleInSLC

%         if ( slcLocationIncrease )
%             if ( prevBestHBSlcLocation >= meanHBSlcLocation(b) )
%                 continue;
%             end
%         else
%             if ( prevBestHBSlcLocation <= meanHBSlcLocation(b) )
%                 continue;
%             end
%         end

        indForCurrHB = indTrigger{b};
        binDataForCurrHB = binData(:,:,indForCurrHB);
        triggerTimeForCurrHB = triggerTimeSLC(indForCurrHB);

        currHBImAtRequiredCardiacPhases = performCineInterpolatiorTriggerTime(binDataForCurrHB, triggerTimeForCurrHB, cardiacPhases);

        % perform the moco
        header = CreateFtkHeaderInfo(currHBImAtRequiredCardiacPhases, [1.8 1.8 3.6]);
        deformInitial = zeros(size(currHBImAtRequiredCardiacPhases));
        
        [moco, dx, dy, dxInv, dyInc] = Matlab_PerformTemporalPairWiseMotionCorrection(double(prevHBImAtRequiredCardiacPhases), double(currHBImAtRequiredCardiacPhases), ...
            header, inverse, initial, iters, sigma, neighbor, stepDiv, moreIterInv, algo, deformInitial, deformInitial, deformInitial, deformInitial, volumePreserving);
        
        dxSqr = dx.*dx;
        dySqr = dy.*dy;
        sumDeform = dxSqr(indUsed(:)) + dySqr(indUsed(:));
        meanTotalDeform(b) = norm(sumDeform);
    end

    bestCardiacPhase = -1;
    minTotalDeform = inf;
    for b=1:numOfCardiacCycleInSLC
        if ( meanTotalDeform(b) ~= inf )
            if ( meanTotalDeform(b) < minTotalDeform )
                bestCardiacPhase = b;
                minTotalDeform = meanTotalDeform(b);
            end
        end
    end    
    
    [minTotalDeform, bestCardiacPhase] = min(meanTotalDeform);

    costFunc = meanTotalDeform;
end
