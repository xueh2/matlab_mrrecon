

function [binAveAll_BreathingGating_MocoAve, ...
        binAveAll_BreathingGating_Linear, ...
        header_binAveAll_BreathingGating, ...
        binAveAll_BreathingGating_ImagePositionPatient, ...
        binAveAll_BreathingGating_ImageOrientationPatient] = ...
            PerformBreathingRetrogatingAfterSliceRecon(reconedSliceLocations, Nfe, Npe, NPhs, binAveAll, sliceThickness, ... 
                        triggerTime, sliceLocation, cardiacPhases, binOverlapRatio, temporalRes, dataAcquired, resporitoryAcceptanceWindow, perfromMTDCombination, ...
                        strategy, inverse, initial, numOfPre, iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving, ...
                        bestHeartBeatsImagePositionPatient, bestHeartBeatsImageOrientationPatient, mocoDir)
% this function perform the breathing retro gating after the slice recon
% both linear interpolation and moco+ave is used to estimate the missing slices
% the resulting volume is recored for patient position and orientation for every slice
% if the parallel sweeping sampling is used, the volume can be directly used
% if the rotation sweeping sampling is used, the volume needs a scatter interpolation step
%
% Inputs:
%    dataDir : folder storing the 4D datasets in dicom format (*.ima or *.dcm files)
%    plotFlag : if 1, then plot the information related to the datasets
%
% Output:
%    data4D: it is a structure including following fields:
%    all fields are stored after the sorting    
%
%    filenames : dicom file name
%    info : dicom info loaded
%    sliceLocation
%    triggerTime
%    acquisitionTime 
%    dataAcquired : loaded image data before correction
%    imagePositionPatient : position patient before correction
%    imageOrientationPatient : orientation patient before correction
%    imagePositionPatientFixed  : position patient after correction
%    imageOrientationPatientFixed  : orientation patient after correction
%    dataAcquiredFixed : image data after correction
%
%     ***************************************
%     *  Hui Xue (hui-xue@siemens.com       *
%     *  2012-04                            *
%     ***************************************

[startSLC, minInd] = min(reconedSliceLocations);
[endSLC, maxInd] = max(reconedSliceLocations);

slices = [startSLC:(sign(endSLC-startSLC)*sliceThickness):endSLC];
NSlc = numel(slices);

binAveAll_BreathingGating_ImagePositionPatient = zeros(NSlc, 3);
binAveAll_BreathingGating_ImageOrientationPatient = zeros(NSlc, 6);

binAveAll_BreathingGating_Linear = zeros(Nfe, Npe, NSlc, NPhs);
binAveAll_BreathingGating_MocoAve = zeros(Nfe, Npe, NSlc, NPhs);
for s=1:NSlc
    
    if ( slices(s) <= startSLC )
        binAveAll_BreathingGating_Linear(:,:,s,:) = binAveAll(:,:,minInd,:);
        binAveAll_BreathingGating_MocoAve(:,:,s,:) = binAveAll(:,:,minInd,:);        
        binAveAll_BreathingGating_ImagePositionPatient(s, :) = bestHeartBeatsImagePositionPatient(minInd, :);
        binAveAll_BreathingGating_ImageOrientationPatient(s, :) = bestHeartBeatsImageOrientationPatient(minInd, :);
        continue;
    end
    
    if ( slices(s) >= endSLC )
        binAveAll_BreathingGating_Linear(:,:,s,:) = binAveAll(:,:,maxInd,:);
        binAveAll_BreathingGating_MocoAve(:,:,s,:) = binAveAll(:,:,maxInd,:);
        binAveAll_BreathingGating_ImagePositionPatient(s, :) = bestHeartBeatsImagePositionPatient(maxInd, :);
        binAveAll_BreathingGating_ImageOrientationPatient(s, :) = bestHeartBeatsImageOrientationPatient(maxInd, :);        
    end
    
    for k=1:numel(reconedSliceLocations)-1
        if ( slices(s)>=reconedSliceLocations(k) & slices(s)<reconedSliceLocations(k+1) )
            break;
        end
    end
    
    alpha = (slices(s)-reconedSliceLocations(k))/(reconedSliceLocations(k+1)-reconedSliceLocations(k));
    binAveAll_BreathingGating_Linear(:,:,s,:) = (1-alpha)*binAveAll(:,:,k,:) + alpha*binAveAll(:,:,k+1,:);
    
    % compute the new image position and orientation
    binAveAll_BreathingGating_ImagePositionPatient(s, :) = (1-alpha)*bestHeartBeatsImagePositionPatient(k,:) + alpha*bestHeartBeatsImagePositionPatient(k+1,:);
    binAveAll_BreathingGating_ImageOrientationPatient(s, 1:3) = (1-alpha)*bestHeartBeatsImageOrientationPatient(k, 1:3) + alpha*bestHeartBeatsImageOrientationPatient(k+1, 1:3);
    binAveAll_BreathingGating_ImageOrientationPatient(s, 4:6) = (1-alpha)*bestHeartBeatsImageOrientationPatient(k, 4:6) + alpha*bestHeartBeatsImageOrientationPatient(k+1, 4:6);
    
    keyFrames = binAveAll_BreathingGating_Linear(:,:,s,:);
    
    % perform another round of retro-moco and averaging
    currSliceRange = [slices(s)-sliceThickness slices(s)+sliceThickness];
    for b=1:NPhs
               
%         indPhs = find( ( triggerTime>=(cardiacPhases(b)-binOverlapRatio*temporalRes) ) ... 
%             & ( triggerTime<(cardiacPhases(b)+binOverlapRatio*temporalRes) ) ...
%             & (sliceLocation>=currSliceRange(1)) & (sliceLocation<currSliceRange(2)) );
%         
%         numOfIm = numel(indPhs);       
%         
%         disp(['b ' num2str(b) ' s ' num2str(s) ' numofIm ' num2str(numOfIm)]);
        
        numOfIm = 0;
        binOverlapRatioUsed = binOverlapRatio;
        
        indPhs = find( ( triggerTime>=(cardiacPhases(b)-binOverlapRatioUsed*temporalRes) ) ... 
                & ( triggerTime<(cardiacPhases(b)+binOverlapRatioUsed*temporalRes) ) ...
                & (sliceLocation>=currSliceRange(1)) & (sliceLocation<currSliceRange(2)) );
        numOfIm = numel(indPhs);
        
        while ( numOfIm < 2 )
        
            indPhs = find( ( triggerTime>=(cardiacPhases(b)-binOverlapRatioUsed*temporalRes) ) ... 
                & ( triggerTime<(cardiacPhases(b)+binOverlapRatioUsed*temporalRes) ) ...
                & (sliceLocation>=currSliceRange(1)) & (sliceLocation<currSliceRange(2)) );

            numOfIm = numel(indPhs);       
        
            binOverlapRatioUsed = 1.1*binOverlapRatioUsed;
        end
        
        if ( numOfIm > 0 )

            binData = zeros([size(dataAcquired(:,:,1,1)) numOfIm+1]);    
            for ii=1:numOfIm
                dataX = dataAcquired(:,:,indPhs(ii));
                binData(:,:,ii) = dataX;
            end            
            
            % get the keyframe image
            currKeyFrame = keyFrames(:,:,b);
            binData(:,:,end) = currKeyFrame;
            
            % perform the moco
            header = CreateFtkHeaderInfo(binData, [1 1 1]);
            % header.sizeZ = size(binData, 3);

            [moco, dx, dy, dxInv, dyInv] = Matlab_PerformTemporalMotionCorrection(binData, header, numOfIm, strategy, inverse, initial, numOfPre, ...
                iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving);

            % perform the average
            % binAve = mean(moco(:,:,1:end-1), 3);
            % perform the average with selected frames
            mocoQuality = zeros(numOfIm, 1);
            header2D = header;
            header2D.sizeZ = 1;            
            for kk=1:numOfIm                
                [meanNorm, maxNorm] = analyzeDeformationField2DNoJacobian(dx(:,:,kk), dy(:,:,kk), header2D);
                mocoQuality(kk) = meanNorm;
            end            
            [mocoQualitySorted, indMoCOQuality] = sort(mocoQuality);            
            numOfImageUsed = ceil(resporitoryAcceptanceWindow*numOfIm);
            if ( numOfImageUsed < 1 ) 
                numOfImageUsed = 1;
            end
            
            if ( numOfImageUsed < 4 )
                numOfImageUsed = 4;
            end
            
            if ( numOfIm <= 4 )
                numOfImageUsed = numOfIm;
            end
            
            indImageUsed = indMoCOQuality(1:numOfImageUsed);
            mocoUsedForAve = moco(:,:,indImageUsed(:));
            dxUsedForAve = dx(:, :, indImageUsed(:));
            dyUsedForAve = dy(:, :, indImageUsed(:));            
            header.sizeZ = numOfImageUsed;
            
            if ( perfromMTDCombination & numOfImageUsed>=4  )
                % the deformation fields weighted combination
                deltaT = 0.05;
                maxIter = 30;
                beta = 0.1;
                mu = 1;
                tol = 1e-4;
                BSplineDeriv = 0;
                [weights, finalIter] = estimateAveragingWeights(mocoUsedForAve, header, dxUsedForAve, dyUsedForAve, deltaT, maxIter, mu, beta, tol, BSplineDeriv);

                ind = find(weights<0);
                weights(ind(:)) = 0;

                weightsSum = sum(weights, 3);
                weightsSum(find(weightsSum<0.01)) = 1.0;

                weights2 = weights;
                for pp=1:header.sizeZ
                    weights2(:,:, pp) = weights(:,:, pp) ./ weightsSum;
                end

                mocoUsedForAve = mocoUsedForAve .* weights2;
                binAve = sum(mocoUsedForAve, 3);                
            else
                binAve = mean(mocoUsedForAve, 3);
            end            
                        
            binAveAll_BreathingGating_MocoAve(:,:,s,b) = binAve;
            
            e = norm(binAve(:));
            if ( e < 1e2 )
                warning('Binned image may be empty ... ');
                figure; imagescn(binAve);
            end
        else
            warning('No image found for this bin ... ');
        end
    end
end

% switch the bin and slice

% compute the image position and orientation for the recon volume
firstSlcImagePosition = bestHeartBeatsImagePositionPatient(minInd, :);
firstSlcImageOrientation = bestHeartBeatsImageOrientationPatient(minInd, :);
lastSlcImagePosition = bestHeartBeatsImagePositionPatient(maxInd, :);

header_binAveAll_BreathingGating = header;
header_binAveAll_BreathingGating.sizeX = size(binAveAll_BreathingGating_MocoAve, 2);
header_binAveAll_BreathingGating.sizeY = size(binAveAll_BreathingGating_MocoAve, 1);
header_binAveAll_BreathingGating.sizeZ = size(binAveAll_BreathingGating_MocoAve, 3);
header_binAveAll_BreathingGating.sizeT = size(binAveAll_BreathingGating_MocoAve, 4);
header_binAveAll_BreathingGating.spacingZ = sliceThickness;
header_binAveAll_BreathingGating.spacingT = temporalRes;
header_binAveAll_BreathingGating.positionPatient = firstSlcImagePosition;
header_binAveAll_BreathingGating.orientationPatient(:,1) = firstSlcImageOrientation(1:3); % row vector, x
header_binAveAll_BreathingGating.orientationPatient(:,2) = firstSlcImageOrientation(4:6); % col vector, y
% norm, z
header_binAveAll_BreathingGating.orientationPatient(:,3) = (lastSlcImagePosition - firstSlcImagePosition) / norm(lastSlcImagePosition - firstSlcImagePosition); % row vector

filename = fullfile(mocoDir, ['binAveAll_BreathingGating_MocoAve']);
Matlab_SaveAnalyze(single(binAveAll_BreathingGating_MocoAve), header_binAveAll_BreathingGating, [filename '.hdr']);

filename = fullfile(mocoDir, ['binAveAll_BreathingGating_Linear']);
Matlab_SaveAnalyze(single(binAveAll_BreathingGating_Linear), header_binAveAll_BreathingGating, [filename '.hdr']);

save Im4D_BreathingGating_Linear binAveAll_BreathingGating_MocoAve binAveAll_BreathingGating_Linear header_binAveAll_BreathingGating
