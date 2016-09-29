


%function [Im4D, cardiacPhases, slices, navigators, others] = PerformCardiac4DReconFloatingSliceLocation_20120307(dataDir, sliceThickness, temporalRes, option)
% function [Im4D, cardiacPhases, slices, navigators, others] = PerformCardiac4DReconFloatingSliceLocation_20120307(dataDir, sliceThickness, temporalRes)
%
% this function performs the slice-by-slice 4D recon for floaing slice
% acqusition datasets.
%
% Inputs:
%    dataDir : folder storing the 4D datasets in dicom format (*.ima or *.dcm files)
%    sliceThickness : required slice thickness in mm for the recon results
%    temporalRes : required temporal resolution in ms for the recon results
%    option.binWidthSLC : the width of bin along SLC, in mm; for a bin center location b, the covered range is [b-binWidthSLC b+binWidthSLC]
%    option.binOverlapRatioTemporal : for a cardiac phase p, the covered range is [p-binOverlapRatioTemporal*temporalRes p+binOverlapRatioTemporal*temporalRes]
%    option.strategyForBestHeartBeat: strategy to find the next best heart beat
%    option.plotFlag : if 1, then plot the information related to the datasets
%
% Output:
%    Im4D : [COL LIN PHS SLC], the 4D recon images
%    cardiacPhases : [PHS 1] vector, storing the cardiac phases reconstructed
%    slices : [SLC 1] vector, storing the slice locations reconstructed
%    navigators : estimated navigators in each SLC
%    others.frames : picked frames for every SLC
%    others. ...

%     ***************************************
%     *  Hui Xue (hui-xue@siemens.com       *
%     *  2012-03                            *
%     ***************************************

%% parameters

% gif parameters
sizeRatio = option.sizeRatio;
centre = option.centre;
width = option.width;
delay = option.delay;

% moco parameters
strategy = 'Direct';
inverse = 1;
initial = 0;
numOfPre = 0;
iters = [64 32 16];
sigma = 12;
neighbor = 2.0;
stepDiv = 3.0;
moreIterInv = 0;
algo = 'GLCC';
volumePreserving = 0;

% binnign parameters
binWidthSLC = option.binWidthSLC;
binOverlapRatio = option.binOverlapRatioTemporal;
plotFlag = option.plotFlag;

%% get the folder information
cd(dataDir)
currDir = dataDir;
mocoDir = fullfile(currDir, option.mocoDir);
mkdir(mocoDir);

[names, numFile] = findFILE(fullfile(currDir, 'dicom'), '*.ima');
if ( numFile == 0 )
    [names, numFile] = findFILE(fullfile(currDir, 'dicom'), '*.dcm');
end

if ( numFile == 0 )
    error(['No images are found at ' currDir '/dicom']);
end

disp('-------------------------------------------------------');

filename = fullfile(currDir, 'dicom', 'data4D.mat');
if ( ~isFileExist(filename) )
    
    data4D = LoadAndCorrectDicomImages(fullfile(currDir, 'dicom'), plotFlag);
    sliceLocation = data4D.sliceLocation;
    triggerTime = data4D.triggerTime;
    acquisitionTime = data4D.acquisitionTime;
    dataAcquired = data4D.dataAcquiredFixed;
    imagePositionPatient = data4D.imagePositionPatientFixed;
    imageOrientationPatient = data4D.imageOrientationPatientFixed;
    
    filenames = cell(numel(data4D),1);
    info = cell(numel(data4D),1);
    
    for k=1:numel(data4D)
        filenames{k} = data4D(k).filenames;
        info{k} = data4D(k).info;
    end
    
    save data4D triggerTime acquisitionTime dataAcquired imagePositionPatient imageOrientationPatient sliceLocation filenames info

else
    load(filename);
end

%% get pixel spacing
info = dicominfo(names{1});
data = dicomread(names{1});
header = CreateFtkHeaderInfo(data, [info.PixelSpacing' info.SliceThickness]);
Nfe = header.sizeX;
Npe = header.sizeY;
pixelSpacing = [info.PixelSpacing' info.SliceThickness];

%% if necessary, load the alternative data and replace the current one
% if ( exist('alterDataName') )
%     if ( ~isempty(alterDataName) )
%         load(alterDataName);
%         dataAcquired = X_display;
%         headerAll = header;
%         headerAll.sizeZ = size(dataAcquired, 3);
%         dataAcquired = dataAcquired*1e5;
%         Matlab_SaveAnalyze(single(dataAcquired), headerAll, ['alterDataAll' '_' option.mocoDir '.hdr']);        
%     end
% end

%% if the slice location is not monotolic, it must be th rotating cases
sliceLocationNew = recomputeSliceLocation(dataAcquired, pixelSpacing, sliceLocation, imagePositionPatient, imageOrientationPatient);

sliceLocationOri = sliceLocation;
sliceLocation = sliceLocationNew;

% %% load the dicom images
% disp('Loading the images ... ');
% 
 
% imagePositionPatient = zeros(numFile, 3);
% imageOrientationPatient = zeros(numFile, 6);
% sliceLocation = zeros(numFile, 1); % in the unit of mm
% triggerTime = zeros(numFile, 1); % in the unit of ms
% acquisitionTime = zeros(numFile, 1);  % in the unit of ms
% dataAcquired = zeros(header.sizeY, header.sizeX, numFile);
% infoes = cell(numFile, 1);
% imageNumber = zeros(numFile, 1);
% filenames = cell(numFile, 1);
% 
% tic
% 
% % load image info
% for k=1:numFile
%     [pathstr, name, ext] = fileparts(names{k});
%     infoes{k} = dicominfo(names{k});
%     imageNumber(k) = infoes{k}.InstanceNumber;
% end
% 
% % sort by imagenumber
% [imageNumber2, ind] = sort(imageNumber);
% 
% % load image content
% for k=1:numFile    
%     [pathstr, name, ext] = fileparts(names{ind(k)});
%     info = infoes{ind(k)};
%     sliceLocation(k) = info.SliceLocation;
%     triggerTime(k) = info.TriggerTime;
%     timeInSeconds = ConvertDicomAcquisitionTime2Seconds(info.AcquisitionTime);
%     acquisitionTime(k) = timeInSeconds;
%     imagePositionPatient(k, :) = info.ImagePositionPatient;
%     imageOrientationPatient(k, :) = info.ImageOrientationPatient;
%     filenames{k} = names{ind(k)};
%     
%     disp([num2str(ind(k)) ' - ' name ' - ' num2str(sliceLocation(k)) ' - ' num2str(triggerTime(k)) ' - ' info.AcquisitionTime ' - ' num2str(timeInSeconds)]);
%     
%     dataX = dicomread(names{ind(k)});    
%     dataAcquired(:,:,k) = dataX;
% end
% 
% headerAll = CreateFtkHeaderInfo(dataAcquired, [info.PixelSpacing' info.SliceThickness]);
% Matlab_SaveAnalyze(single(dataAcquired), headerAll, fullfile(currDir, 'dataAcquired.hdr'));
% disp(['Loading images : ' num2str(toc)]);
% 
% %% check and fix the data properties
% 
% % check to fix the image orientation
% pixelSpacing = [info.PixelSpacing' info.SliceThickness];
% [imagePositionPatientFixed, imageOrientationPatientFixed, dataAcquiredFixed] = fixImageOrientation(imagePositionPatient, imageOrientationPatient, dataAcquired, pixelSpacing);
% 
% save dataReadyForProcessing filenames info sliceLocation triggerTime acquisitionTime imagePositionPatient ... 
%     imageOrientationPatient imagePositionPatientFixed imageOrientationPatientFixed dataAcquiredFixed dataAcquired headerAll
% 
% headerAll = CreateFtkHeaderInfo(dataAcquiredFixed, [info.PixelSpacing' info.SliceThickness]);
% Matlab_SaveAnalyze(single(dataAcquiredFixed), headerAll, fullfile(currDir, 'dataAcquiredFixed.hdr'));
% 
% % if necessary, correct acquisition time
% acquisitionTime = acquisitionTime - acquisitionTime(1);
% acquisitionDelta = acquisitionTime(2:end) - acquisitionTime(1:end-1);
% ind = find(acquisitionDelta>1)
% if ( ~isempty(ind) )
%     
%     numOfGaps = numel(ind);
%     acquisitionTime2 = acquisitionTime;
%     
%     for g=1:numOfGaps
%         gap = acquisitionTime(ind(g)+1) - acquisitionTime(ind(g))
% 
%         gap2 = acquisitionTime(ind(g)) - acquisitionTime(ind(g)-1)
% 
%         acquisitionTime2(ind(g)+1:end) = acquisitionTime2(ind(g)+1:end) - gap + gap2;
%     end
%     
%     figure;
%     scatter(1:numFile, acquisitionTime2, '+');
%     title([currDir '- ' num2str(numFile) ' images']);
%     box on
%     xlabel('Image number');
%     ylabel('Acquisition time corrected');
%     
%     acquisitionTimeOri = acquisitionTime;
%     acquisitionTime = acquisitionTime2;
% end
% 
if ( plotFlag )
    figure;
    scatter(triggerTime, sliceLocation, '.');
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Trigger Time - ms');
    ylabel('Slice Location - mm');

    figure;
    scatter(imageNumber2, triggerTime, '.');
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Image number');
    ylabel('Trigger time - ms');

    figure;
    scatter(imageNumber2, sliceLocation, '.');
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Image number');
    ylabel('Slice location - mm');

    figure;
    scatter(imageNumber2, acquisitionTime, '.');
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Image number');
    ylabel('Acquisition time');

    figure;
    scatter(sliceLocation, acquisitionTime, '.');
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Slice Location - mm');
    ylabel('Acquisition Time - ms');
    
    figure;
    hold on
    scatter3(imagePositionPatient(:, 1), imagePositionPatient(:,2), imagePositionPatient(:,3));
    hold off
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('ImagePositionPatientX - mm');
    ylabel('ImagePositionPatientY - mm');
    zlabel('ImagePositionPatientZ - mm');
    
    figure;
    hold on
    scatter3(imageOrientationPatient(:, 1), imageOrientationPatient(:,2), imageOrientationPatient(:,3));
    hold off
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('imageOrientationPatient row vector x - mm');
    ylabel('imageOrientationPatient row vector y - mm');
    zlabel('imageOrientationPatient row vector z - mm');
        
    plotStep = 30;
    ha = -1;
    %aviobj = avifile('slice3D.avi');
    for k=1:plotStep:numFile
        [d, h] = LoadDicomImageMrFtk(filenames{k});
        ha = plotMrFtkImageIn3D(d, h, ha, 200, 400);
        frame = getframe(gcf);
        %aviobj = addframe(aviobj,frame);
    end
    %aviobj = close(aviobj);
end

%% compute the mean R-R interval, mean starting phase and mean ending phase and other cardiac phases
imageNumber2 = 1:numel(sliceLocation);

disp('-----------------------------------------');
% all indexes are in the order of global frames
[meanRR, meanTriggerTimeDelta, startingPhase, startingInd, endingPhase, endingInd, indexes] = ComputeMeanRR(triggerTime);
if ( plotFlag )
    figure;
    hold on
    scatter(1:numel(imageNumber2), triggerTime, '.');
    plot(startingInd, startingPhase, 'rx');
    plot(endingInd, endingPhase, 'gx');    
    hold off
    title([currDir '- ' num2str(numFile) ' images']);
    box on
    xlabel('Image number');
    ylabel('Trigger time - ms');
end

Nfe = size(dataAcquired, 1);
Npe = size(dataAcquired, 2);

numOfCompleteRR = numel(startingPhase);

startPhase = mean(startingPhase);
endPhase = mean(endingPhase);

if ( ~exist('cardiacPhases') )
    cardiacPhases = [startPhase:temporalRes:endPhase];
else
    disp(['Preset cardiac phases ... ']);
end

NPhs = numel(cardiacPhases);
headerPhs = header;
headerPhs.sizeZ = NPhs;

disp(['Mean RR : ' num2str(meanRR)]);
disp(['Mean trigger time delta : ' num2str(meanTriggerTimeDelta)]);
disp(['Number of complete cardiac cycle : ' num2str(numOfCompleteRR)]);
disp(['Number of reconstructed cardiac phases : ' num2str(numel(cardiacPhases))]);

%% compute the mean slice location delta
sliceLocationDelta = sliceLocation(2:end) - sliceLocation(1:end-1);
meanSliceLocationDelta = mean(abs(sliceLocationDelta));
disp(['Mean slice location delta in acqusition : ' num2str(meanSliceLocationDelta)]);

if ( sliceThickness < 4*meanSliceLocationDelta )
    warning(['Required slice thickenss ' num2str(sliceThickness) ' is too small, compared to floating stepsize ' num2str(meanSliceLocationDelta) ' in this dataset']);
    sliceThickness = 4*meanSliceLocationDelta;
end

% find the slice location is going to increase or decrease
if ( mean(sliceLocationDelta) == 0 )
    error('Slice location does not change between slices ... ');
end

slcLocationIncrease = 1;
if ( mean(slcLocationIncrease) < 0 )
    slcLocationIncrease = 0;
end

startSLC = min(sliceLocation);
startSLC = sign(startSLC) * (abs(startSLC)+sliceThickness/2);

endSLC = max(sliceLocation);
endSLC = sign(endSLC) * (abs(endSLC)-sliceThickness/2);

slices = [startSLC:(sign(endSLC-startSLC)*sliceThickness):endSLC];
NSlc = numel(slices);
disp(['Number of reconstructed slices : ' num2str(numel(slices))]);

% check the slice window width
acquisitionDelta = acquisitionTime(2:end) - acquisitionTime(1:end-1);
meanSlcAcqDuration = mean(acquisitionDelta);

SlcAcqTimeInSlidingWindow = 2*meanSlcAcqDuration*binWidthSLC/meanSliceLocationDelta;
disp(['Approx. slice acquisition time in sliding window : ' num2str(SlcAcqTimeInSlidingWindow)]);
disp(['Approx. breathing cycles within the sliding window : ' num2str(SlcAcqTimeInSlidingWindow/3)]);

if ( SlcAcqTimeInSlidingWindow/3 < 2 )
    warning('The slice window width is too narrow ...');
    binWidthSLC = 3*meanSliceLocationDelta/meanSlcAcqDuration;
    disp(['Adjust slice sliding window width : ' num2str(binWidthSLC)]);
    SlcAcqTimeInSlidingWindow = 2*meanSlcAcqDuration*binWidthSLC/meanSliceLocationDelta;
end

disp(['Approx. number of heart beats within the SLC sliding window : ' num2str(SlcAcqTimeInSlidingWindow*1000/meanRR)]);

%% perform the slice by slice recon
disp('-----------------------------------------');

% strategy three : dijkstra's algorithm
binAveAll = zeros(Nfe, Npe, NSlc, NPhs);
suffix = ['_' option.strategyForBestHeartBeat];
navigatorAll = cell(NSlc, 1);

% starting from the middle SLC
midSLC = floor(median(1:NSlc));

others.frames = cell(NSlc, 1);

prevBestHB = 0; % HB : heart beat
prevBestHBIm = [];
prevBestHBSlcLocation = [];
mask = [];
if ( isFileExist('mask.mat') )
    load('mask.mat');
end

slcInd = midSLC:NSlc;
slcInd = [slcInd midSLC:-1:1];

crossHalfSLCFlag = 1;

bestHeartBeats = zeros(NSlc, 1);
reconedSliceLocations = zeros(NSlc, 1);
bestHeartBeatsKeyFrame = zeros(NSlc, 1);
bestHeartBeatsImagePositionPatient = zeros(NSlc, 3);
bestHeartBeatsImageOrientationPatient = zeros(NSlc, 6);

numOfImageInEachBins = zeros(NPhs, NSlc);

for currSlc=1:NSlc+1
    
    s = slcInd(currSlc);
    
    disp(['Recon slice : ' num2str(s)]);
    
    if ( s==midSLC )
        currSliceLocation = slices(s);
    else
        currSliceLocation = prevBestHBSlcLocation;
    end
    
    sliceRange = [currSliceLocation-binWidthSLC currSliceLocation+binWidthSLC];
    
    % find the images falling into this SLC
    ind = find( sliceLocation>=sliceRange(1) & (sliceLocation<sliceRange(2)) );
    
    % find the complete heart beat that this SLC involved
    % ind is modified to match the entire heart beats
    [indHB, ind] = findCompleteHeartBeats(ind, indexes);
    
    % if the recon has reached the end
    if ( numel(indHB) < 4 )
        continue;
    end
    
    % get the bin data
    numOfIm = numel(ind)
    binData = zeros([Nfe Npe numOfIm]);    
    for ii=1:numOfIm
        dataX = dataAcquired(:,:,ind(ii));
        binData(:,:,ii) = dataX;
    end

    others.frames{s} = ind;
    
    if ( isempty(mask) )
        h = figure; imshow(mean(binData, 3), [], 'InitialMagnification', 500);
        BW = roipoly;
        close(h);
        mask = BW;
        save('mask.mat', 'mask');
    end
    
    % find the total heart beats in this SLC
    triggerTimeSLC = triggerTime(ind);
    % [RR, TD, sPhs, sInd, ePhs, eInd, indTrigger] = ComputeMeanRR(triggerTime(ind));
    [RR, TD, sPhs, sInd, ePhs, eInd, indTrigger] = ComputeMeanRRForCompleteHB(indHB, indexes, triggerTime);
    if ( plotFlag )
        figure;
        hold on
        scatter(sInd(1):eInd(end), triggerTime(ind), '.');
        plot(sInd, sPhs, 'rx');
        plot(eInd, ePhs, 'gx');    
        hold off
        title([num2str(numel(sPhs)) ' heart beats']);
        box on
        xlabel('Image number');
        ylabel('Trigger time - ms');
    end

    numOfCardiacCycleInSLC = numel(sPhs);
    
    acquisitionTimeSLC = acquisitionTime(ind);
    
    % compute the image indexes for every cardiac cycle in the original image orders
    indTriggerInOriginalOrder = indTrigger;
%     for b=1:numOfCardiacCycleInSLC        
%         currInd = indTrigger{b};
%         currIndOri = indTriggerInOriginalOrder{b};
%         for ii=1:numel(currInd)
%             currIndOri(ii) = ind(currInd(ii));
%         end
%         indTriggerInOriginalOrder{b} = currIndOri;
%     end
    
    % compute the mean slice location for every heart beats
    meanHBSlcLocation = zeros(numOfCardiacCycleInSLC, 1);
    for b=1:numOfCardiacCycleInSLC 
        slcLoc = sliceLocation( indTriggerInOriginalOrder{b} );
        meanHBSlcLocation(b) = mean(slcLoc);
    end
    
    %% detect the best heart beats given a SLC    
    if ( s==midSLC )   
        % for the middle SLC, use the navigator signal
        
        crossHalfSLCFlag = ~crossHalfSLCFlag;        
        if ( crossHalfSLCFlag )
            slcLocationIncrease = ~slcLocationIncrease;
        end
        
        % find the navigator signal
        headerNav = header;
        headerNav.sizeZ = numOfIm;
        navigator = computeImageBasedNavigator(binData, headerNav, acquisitionTime(ind), triggerTime(ind), mask, 24, [64 64 64]);
        navigators{s} = navigator;      

        % bin along triggerTimes
        triggerTimeSLC = triggerTime(ind);    
        range = max(triggerTimeSLC) - min(triggerTimeSLC)
        binNumber = floor(0.9*range/temporalRes)
        [N, X] = hist(triggerTimeSLC, binNumber);

        % find the best heart beat    
        totalImageInBins = 0;
        imageInBins = cell(numel(X)-1, 1);
        for b=1:numel(X)-1

            indSlcTime = find( (triggerTime>=(X(b)-temporalRes))&(triggerTime<(X(b)+temporalRes)) ...
                & (sliceLocation>=sliceRange(1))&(sliceLocation<sliceRange(2))  );

            imageInBins{b} = indSlcTime;

            if ( isempty(indSlcTime) )
                imageInBins{b} = imageInBins{b-1};
            end

            totalImageInBins = totalImageInBins + numel(indSlcTime);
        end

        % -----------------------------------
        % find the key frames
        bestKeyFrameForEachBin = zeros(numel(X)-1, 1);

        costAndPredecessor = cell(numel(X)-1, 1); % every element is a [numOfImageInBins 2] array, storing (predessor, bestCost)
        costAndPredecessor{1} = zeros(numel(imageInBins{1}), 2);

        for b=2:numel(X)-1

            imPreInd = imageInBins{b-1};
            imInd = imageInBins{b};

            costAndPredecessor{b} = zeros(numel(imInd), 2);

            preBestCost = costAndPredecessor{b-1};
            bestCost = zeros(numel(imInd), 2);
            for im=1:numel(imInd)

                indNav = find(ind==imInd(im));
                nav = navigator(indNav);

                minCost = inf;
                for imPre=1:numel(imPreInd)

                    indPreNav = find(ind==imPreInd(imPre));
                    navPre = navigator(indPreNav);

                    costForEdge = abs(navPre-nav) + preBestCost(imPre, 2);

                    if ( costForEdge < minCost )
                        minCost = costForEdge;
                        bestCost(im, 1) = imPreInd(imPre);
                        bestCost(im ,2) = minCost;
                    end
                end
            end

            costAndPredecessor{b} = bestCost;
        end

        for b=numel(X)-1:-1:1

            imInd = imageInBins{b};
            bestCost = costAndPredecessor{b};

            if ( b == numel(X)-1 )        

                preInd = bestCost(:, 1);
                minCost = bestCost(:, 2);

                [v, keyFrame] = min(minCost);
                bestKeyFrameForEachBin(b) = imInd(keyFrame);
                preInd = preInd(keyFrame);
            else
                bestKeyFrameForEachBin(b) = preInd;        
                pInd = find(imInd==preInd);
                preInd = bestCost(pInd, 1);
            end
        end

        % find the best heart beat as the heart beat containing the most number of selected keyframes
        voteForKeyFrame = zeros(numOfCardiacCycleInSLC, 1);    
        for b=1:numel(X)-1
            keyFrame = bestKeyFrameForEachBin(b);
            for v=1:numOfCardiacCycleInSLC
                indCurrPhase = indTriggerInOriginalOrder{v};
                indKeyFrame = find(keyFrame==indCurrPhase);
                if ( ~isempty(indKeyFrame) )
                    voteForKeyFrame(v) = voteForKeyFrame(v) + 1;
                end
            end
        end

        [vote, bestCardiacCycle] = max(voteForKeyFrame);
        
        if ( plotFlag )
            figure;
            hold on
            plot(acquisitionTimeSLC, navigator);
            plot(acquisitionTimeSLC, navigator, 'x');
            for bb=1:numOfCardiacCycleInSLC
                plot(acquisitionTimeSLC(sInd(bb)-ind(1)+1), navigator(sInd(bb)-ind(1)+1), 'r+');
                plot(acquisitionTimeSLC(eInd(bb)-ind(1)+1), navigator(eInd(bb)-ind(1)+1), 'r+');
            end
            plot(acquisitionTimeSLC(sInd(bestCardiacCycle)-ind(1)+1), navigator(sInd(bestCardiacCycle)-ind(1)+1), 'g+');
            plot(acquisitionTimeSLC(eInd(bestCardiacCycle)-ind(1)+1), navigator(eInd(bestCardiacCycle)-ind(1)+1), 'g+');
            hold off
            title([' Best heart beat picked - initial']);
            box on
            xlabel('Image number');
            ylabel('Navigator - mm');
        end        
    else
        % for other SLCs, propagate the best cardiac phases from its previous SLC
        
        % compute the interpolated images at the required cardiac phases
        prevHBImAtRequiredCardiacPhases = performCineInterpolatiorTriggerTime(prevBestHBIm, prevBestHBTriggerTime, cardiacPhases, 'linear');
        
        [bestCardiacCycle, costFunc] = findBestHeartBeat(cardiacPhases, meanHBSlcLocation, dataAcquired, indTrigger, triggerTime, ... 
                    prevHBImAtRequiredCardiacPhases, prevBestHBSlcLocation, slcLocationIncrease, option.strategyForBestHeartBeat, mask);
        
%         %% strategy 1 : the mean SSD
%         % compute the mean of total SSD between every current 
%         meanTotalSSD = zeros(numOfCardiacCycleInSLC, 1);
%         meanTotalSSD(:) = -1;
%         
%         % for every current heart beat, first compare its mean slice location to the prev best heart beat
%         % if the slice location moves in the correct direction, this heart beat is lised as candidates
%         % To compute the SSD, first compute the interpolated images of the previous best heart beat at the required cardiac phases
%         % then compute the interpolated images of the current heart beat at the required cardiac phases
%         % the compute the mean SSD between each pair of corresponding cardiac phases
%         
%         % for every current HB
%         for b=1:numOfCardiacCycleInSLC
%             
%             if ( slcLocationIncrease )
%                 if ( prevBestHBSlcLocation >= meanHBSlcLocation(b) )
%                     continue;
%                 end
%             else
%                 if ( prevBestHBSlcLocation <= meanHBSlcLocation(b) )
%                     continue;
%                 end
%             end
%             
%             indForCurrHB = indTrigger{b};
%             binDataForCurrHB = binData(:,:,indForCurrHB);
%             triggerTimeForCurrHB = triggerTimeSLC(indForCurrHB);
%     
%             currHBImAtRequiredCardiacPhases = performCineInterpolatiorTriggerTime(binDataForCurrHB, triggerTimeForCurrHB, cardiacPhases);
%             
%             diffHB = prevHBImAtRequiredCardiacPhases-currHBImAtRequiredCardiacPhases;
%             meanTotalSSD(b) = norm(diffHB(:));
%         end
%         
%         bestCardiacCycle = -1;
%         minTotalSSD = inf;
%         for b=1:numOfCardiacCycleInSLC
%             if ( meanTotalSSD(b) ~= -1 )
%                 if ( meanTotalSSD(b) < minTotalSSD )
%                     bestCardiacCycle = b;
%                     minTotalSSD = meanTotalSSD(b);
%                 end
%             end
%         end
        
        if ( bestCardiacCycle == -1 )
            warning(['SLC ' num2str(s)  ' cannot find a proper heart beat ... ']);
            
            for b=1:numOfCardiacCycleInSLC
                indForCurrHB = indTrigger{b};
                binDataForCurrHB = dataAcquired(:,:,indForCurrHB);
                triggerTimeForCurrHB = triggerTime(indForCurrHB);

                currHBImAtRequiredCardiacPhases = performCineInterpolatiorTriggerTime(binDataForCurrHB, triggerTimeForCurrHB, cardiacPhases, 'linear');

                diffHB = prevHBImAtRequiredCardiacPhases(:)-currHBImAtRequiredCardiacPhases(:);
                meanTotalSSD(b) = norm(diffHB);
            end 
            
            [minTotalSSD, bestCardiacCycle] = min(meanTotalSSD);
        end
    end
        
    % record the reconed slice locations
    reconedSliceLocations(s) = meanHBSlcLocation(bestCardiacCycle);    
    bestHeartBeats(s) = indHB(bestCardiacCycle);
    
    if ( reconedSliceLocations(s) == prevBestHBSlcLocation )
        warning('does not move on any more ...');
    end
    
    % perform the interpolator to compute key frames
    keyFrames = zeros(Nfe, Npe, NPhs);
    
    indForBestCardiacCycle = indTrigger{bestCardiacCycle};
    binDataForBestCardiacCycle = dataAcquired(:,:,indForBestCardiacCycle);
    triggerTimeForBestCardiacCycle = triggerTime(indForBestCardiacCycle);
    
    prevBestHB = bestCardiacCycle;
    prevBestHBIm = binDataForBestCardiacCycle;
    prevBestHBTriggerTime = triggerTimeForBestCardiacCycle;
    prevBestHBSlcLocation = meanHBSlcLocation(bestCardiacCycle);    
    
    currSliceLocation = meanHBSlcLocation(bestCardiacCycle);
    % currSliceRange = [currSliceLocation-binWidthSLC/4 currSliceLocation+binWidthSLC/4];
    currSliceRange = [currSliceLocation-sliceThickness currSliceLocation+sliceThickness];
    
    % correct residual bulk motion in the binData of best cardiac cycle
    header.sizeZ = size(binDataForBestCardiacCycle, 3);
    
    Matlab_SaveAnalyze(single(binDataForBestCardiacCycle), header, fullfile(mocoDir, ['BestBin_SLC' num2str(s) '.hdr']));
    
    if ( option.mocoBestHeartBeat )
        iters_KeyFrame = [32 32 16];
        sigma_KeyFrame = header.sizeY/6+header.sizeX/6;
        [binDataForBestCardiacCycle_moco, dx, dy, invDx, invDy, kf] = PerformMoCoWithKeyFrameSelectionSSD(binDataForBestCardiacCycle, header, strategy, inverse, ...
                                                                                    initial, numOfPre, iters_KeyFrame, sigma_KeyFrame, neighbor, ...
                                                                                    stepDiv, moreIterInv, algo, volumePreserving);    

        Matlab_SaveAnalyze(single(binDataForBestCardiacCycle_moco), header, fullfile(mocoDir, ['BestBin_moco_SLC' num2str(s) '.hdr']));
    else
        kf = PerformKeyFrameSelectionSSD(binDataForBestCardiacCycle);
        binDataForBestCardiacCycle_moco = binDataForBestCardiacCycle;
    end
    
    keyFrames = performCineInterpolatiorTriggerTime(binDataForBestCardiacCycle_moco, triggerTimeForBestCardiacCycle, cardiacPhases, 'pchip');
    Matlab_SaveAnalyze(single(keyFrames), headerPhs, fullfile(mocoDir, ['KeyFrames_SLC' num2str(s) '.hdr']));
    hdr2gifWithWindowSetting(fullfile(mocoDir, ['KeyFrames_SLC' num2str(s) '.hdr']), fullfile(mocoDir, ['KeyFrames_SLC' num2str(s) '.gif']), sizeRatio, centre, width, delay);
    
    % record the image position and orientation for the current best heart beat
    bestHBKeyFrameInd = indexes{bestHeartBeats(s)};
    bestHBKeyFrameInd = bestHBKeyFrameInd(kf+1);
    bestHeartBeatsKeyFrame(s) = bestHBKeyFrameInd;
    bestHeartBeatsImagePositionPatient(s, :) = imagePositionPatient(bestHBKeyFrameInd, :);
    bestHeartBeatsImageOrientationPatient(s, :) = imageOrientationPatient(bestHBKeyFrameInd, :);
    
    % -----------------------------------
    
    for b=1:NPhs
        
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
        
        disp(['b ' num2str(b) ' s ' num2str(s) ' numofIm ' num2str(numOfIm)]);
        
        numOfImageInEachBins(b, s) = numOfIm;
        
        if ( numOfIm == 0 )
            disp('numOfIm == 0');
        end
        
        if ( numOfIm > 0 )

            filename = fullfile(mocoDir, ['bin_triggerTime' num2str(b) '_SliceLocation' num2str(s)]);        
%             [binData, header] = Matlab_LoadAnalyze([filename '.hdr']);

            binData = zeros([size(data) numOfIm+1]);    
            for ii=1:numOfIm
                dataX = dataAcquired(:,:,indPhs(ii));
                binData(:,:,ii) = dataX;
            end            
            
            % get the keyframe image
            currKeyFrame = keyFrames(:,:,b);
            binData(:,:,end) = currKeyFrame;
            
            % perform the moco
            header.sizeZ = size(binData, 3);
%           Matlab_SaveAnalyze(single(binData), header, [filename suffix '.hdr']);

            [moco, dx, dy, dxInv, dyInv] = Matlab_PerformTemporalMotionCorrection(binData, header, numOfIm, strategy, inverse, initial, numOfPre, ...
                iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving);
%             Matlab_SaveAnalyze(single(moco), header, [filename suffix '_moco.hdr']);
%             hdr2gifWithWindowSetting([filename suffix '.hdr'], [filename suffix '.gif'], sizeRatio, centre, width, delay);        

            % -------------------------------------
            % perform the average with selected frames
            mocoQuality = zeros(numOfIm, 1);
            header2D = header;
            header2D.sizeZ = 1;
            
            for kk=1:numOfIm                
                [meanNorm, maxNorm] = analyzeDeformationField2DNoJacobian(dx(:,:,kk), dy(:,:,kk), header2D);
                mocoQuality(kk) = meanNorm;
            end
            
            [mocoQualitySorted, indMoCOQuality] = sort(mocoQuality);
            
            numOfImageUsed = ceil(option.resporitoryAcceptanceWindow*numOfIm);
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
            
            if ( option.useMTDCombination & numOfImageUsed>=4  )
                % ----------------------------------------
                % the deformation fields weighted combination
                deltaT = 0.05;
                maxIter = 30;
                beta = 0.1;
                mu = 1;
                tol = 1e-4;
                BSplineDeriv = 0;
                [weights, finalIter] = estimateAveragingWeights(mocoUsedForAve, header, dxUsedForAve, dyUsedForAve, deltaT, maxIter, mu, beta, tol, BSplineDeriv);
                Matlab_SaveAnalyze(single(weights), header, fullfile(mocoDir,'moco-averaging-weights.hdr'));

                ind = find(weights<0);
                weights(ind(:)) = 0;
                Matlab_SaveAnalyze(single(weights), header, fullfile(mocoDir,'moco-averaging-weights-abs.hdr'));

                weightsSum = sum(weights, 3);
                weightsSum(find(weightsSum<0.01)) = 1.0;

                weights2 = weights;
                for pp=1:header.sizeZ
                    weights2(:,:, pp) = weights(:,:, pp) ./ weightsSum;
                end

                % ------------------------------------------
                % perform the averaging
                mocoUsedForAve = mocoUsedForAve .* weights2;
                binAve = sum(mocoUsedForAve, 3);       
            else
                binAve = mean(mocoUsedForAve, 3);                                     
            end
            
            binAveAll(:,:,s,b) = binAve;
            
            e = norm(binAve(:));
            if ( e<1e3)
                disp('e<1e3');
                figure; imagescn(binAve);
            end            
        else
            disp('No images are found for this bin ... ');
        end
    end
    
    binSlc = binAveAll(:,:,s,:);
    header.sizeZ = size(binSlc, 3);
    filename = fullfile(mocoDir, ['bin_SliceLocation' num2str(s) '_' suffix '_ave']);
    Matlab_SaveAnalyze(single(binSlc), header, [filename '.hdr']);
    hdr2gifWithWindowSetting([filename '.hdr'], [filename '.gif'], sizeRatio, centre, width, delay);
end

binAveAll = single(binAveAll);
save(fullfile(mocoDir, 'HeartBeatPickingResults'), 'bestHeartBeats', ...
    'reconedSliceLocations', 'binAveAll', 'bestHeartBeatsKeyFrame', ...
    'bestHeartBeatsImagePositionPatient', 'bestHeartBeatsImageOrientationPatient', 'numOfImageInEachBins');

header_binAveAll = header;
header_binAveAll.sizeZ = size(binAveAll, 3);
header_binAveAll.sizeT = size(binAveAll, 4);
header_binAveAll.spacingZ = sliceThickness;
header_binAveAll.spacingT = temporalRes;

filename = fullfile(mocoDir, ['binAveAll']);
Matlab_SaveAnalyze(single(binAveAll), header_binAveAll, [filename '.hdr']);

%% ---------------------------------------------------
% one extra moco on the first phase of binAveAll and apply to all other phases, aiming to correct residual mismatch
% strategy = 'Consecutive';
% inverse = 1;
% initial = 0;
% numOfPre = 0;
% iters = [32 16 1];
% sigma = header.sizeY/4;
% neighbor = 2.0;
% stepDiv = 3.0;
% moreIterInv = 1;
% algo = 'GLCC';
% volumePreserving = 0;
% 
% binAveAll_repoCorr = binAveAll;
% 
% % find the best slice as the key frames
% ssd = zeros(NSlc, NSlc);
% for s=1:NSlc    
%     binAveAllSlc = binAveAll(:,:,s,:);    
%     for s2=1:NSlc
%         if ( s <= s2 )
%             binAveAllSlc2 = binAveAll(:,:,s2,:);
%             diffSlcIm = binAveAllSlc-binAveAllSlc2;
%             diffSlcIm = diffSlcIm(:);
%             ssd(s, s2) = sqrt(diffSlcIm'*diffSlcIm);
%             if ( s ~= s2 ) ssd(s2, s) = ssd(s, s2); end
%         end
%     end        
% end
% 
% ssd = sort(ssd, 1);
% mssd = median(ssd);
% [minssd, bestSLC] = min(mssd);
% bestSLC = bestSLC - 1;
% 
% phsForMoCo = 1;
% binAveAllTrigger = binAveAll(:,:,:,phsForMoCo);
% 
% filename = fullfile(mocoDir, ['binAveAll_Phs' num2str(phsForMoCo)]);
% Matlab_SaveAnalyze(single(binAveAllTrigger), headerPhs, [filename '.hdr']);
% 
% headerPhs = header_binAveAll;
% headerPhs.sizeT = 1;
% 
% [moco, dx, dy, dxInv, dyInv] = Matlab_PerformTemporalMotionCorrection(double(binAveAllTrigger), headerPhs, bestSLC, strategy, inverse, initial, numOfPre, ...
%     iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving);
% 
% filename = fullfile(mocoDir, ['binAveAll_MoCo_Phs' num2str(phsForMoCo)]);
% Matlab_SaveAnalyze(single(moco), headerPhs, [filename '.hdr']);
% 
% binAveAll_repoCorr = binAveAll;
% for b=1:NPhs
%     binAveAllTrigger = binAveAll(:,:,:,b);    
%     warpped = Matlab_PerformWarpingSeries2D(double(binAveAllTrigger), headerPhs, single(dx), single(dy), -1, 'BSpline', 5, 1);    
%     binAveAll_repoCorr(:,:,:,b) = warpped;
% end
% 
% filename = fullfile(mocoDir, ['binAveAll_repoCorr']);
% Matlab_SaveAnalyze(single(binAveAll_repoCorr), header_binAveAll, [filename '.hdr']);
% 
% binAveAll_Ori = binAveAll;
% binAveAll = binAveAll_repoCorr;

%% one extra moco and rejection to get rid of mis-matched images

if( 0 )
    size(binAveAll)

    plotStep = 1;
    saveAvi = 1;
    aviName = 'd:/temp/3DRadialVolume.avi';

    slc = 25;
    Render4DVolume(binAveAll, pixelSpacing, bestHeartBeatsImagePositionPatient, bestHeartBeatsImageOrientationPatient, slc, -1, plotStep, 200, 400, 1/20);

    phs = 6;
    Render4DVolume(binAveAll, pixelSpacing, bestHeartBeatsImagePositionPatient, bestHeartBeatsImageOrientationPatient, -1, phs, plotStep, 200, 400, 1/20, saveAvi, aviName);

    figure; plot3(bestHeartBeatsImagePositionPatient(:,1), bestHeartBeatsImagePositionPatient(:,2), bestHeartBeatsImagePositionPatient(:,3), 'x');

    save4DVolume(binAveAll, header_binAveAll_BreathingGating, mocoDir, 'Ori', sizeRatio, centre, width, delay);
end

pixelSpacingRecon = pixelSpacing;
pixelSpacingRecon(3) = sliceThickness;
pixelSpacingRecon

% numOfSubdivision = 8;
% dstPixelSpacing = pixelSpacingRecon;
% 
% [volumeBFFD, headerBFFD] = PerformScatterInterpolationFFD(binAveAll, pixelSpacingRecon, ...
%     bestHeartBeatsImagePositionPatient, bestHeartBeatsImageOrientationPatient, ... 
%     dstPixelSpacing, numOfSubdivision, 'D:/temp/4DRecon/');

%% perform the breathing retrospective gating

strategy = 'Direct';
inverse = 1;
initial = 0;
numOfPre = 0;
iters = [64 32 16];
sigma = 12;
neighbor = 2.0;
stepDiv = 3.0;
moreIterInv = 0;
algo = 'GLCC';
volumePreserving = 0;

[binAveAll_BreathingGating_MocoAve, ...
        binAveAll_BreathingGating_Linear, ...
        header_binAveAll_BreathingGating, ...
        binAveAll_BreathingGating_ImagePositionPatient, ...
        binAveAll_BreathingGating_ImageOrientationPatient] = ...
            PerformBreathingRetrogatingAfterSliceRecon(reconedSliceLocations, Nfe, Npe, NPhs, binAveAll, sliceThickness, ... 
                        triggerTime, sliceLocation, cardiacPhases, binOverlapRatio, temporalRes, dataAcquired, option.resporitoryAcceptanceWindow, option.useMTDCombination, ...
                        strategy, inverse, initial, numOfPre, [100 64 16], 6, neighbor, stepDiv, moreIterInv, algo, volumePreserving, ...
                        bestHeartBeatsImagePositionPatient, bestHeartBeatsImageOrientationPatient, mocoDir);

if 0
    figure; plot3(binAveAll_BreathingGating_ImagePositionPatient(:,1), binAveAll_BreathingGating_ImagePositionPatient(:,2), binAveAll_BreathingGating_ImagePositionPatient(:,3), 'x');

    slc = 35;
    Render4DVolume(binAveAll_BreathingGating_Linear, pixelSpacing, binAveAll_BreathingGating_ImagePositionPatient, binAveAll_BreathingGating_ImageOrientationPatient, slc, -1, 1, 200, 400, 1/20);

    Render4DVolume(binAveAll_BreathingGating_MocoAve, pixelSpacing, binAveAll_BreathingGating_ImagePositionPatient, binAveAll_BreathingGating_ImageOrientationPatient, slc, -1, 1, 200, 400, 1/20);

    phs = 6;
    Render4DVolume(binAveAll_BreathingGating_MocoAve, pixelSpacing, binAveAll_BreathingGating_ImagePositionPatient, binAveAll_BreathingGating_ImageOrientationPatient, -1, phs, 1, 200, 400, 1/20, 0, []);
end

binAveAll_BreathingGating_Linear = single(binAveAll_BreathingGating_Linear);
binAveAll_BreathingGating_MocoAve = single(binAveAll_BreathingGating_MocoAve);

save(fullfile(mocoDir, 'AfterBreathingRetrogatingResults.mat'), 'binAveAll_BreathingGating_Linear', 'binAveAll_BreathingGating_MocoAve', ...
    'header_binAveAll_BreathingGating', 'binAveAll_BreathingGating_ImagePositionPatient', 'binAveAll_BreathingGating_ImageOrientationPatient');
% save AfterBreathingRetrogatingResults binAveAll_BreathingGating_Linear binAveAll_BreathingGating_MocoAve header_binAveAll_BreathingGating binAveAll_BreathingGating_ImagePositionPatient binAveAll_BreathingGating_ImageOrientationPatient

% selectedPhase = 1;
% [binAveAll_BreathingGating_MocoAve, ...
%         binAveAll_BreathingGating_Linear, ...
%         header_binAveAll_BreathingGating, ...
%         binAveAll_BreathingGating_ImagePositionPatient, ...
%         binAveAll_BreathingGating_ImageOrientationPatient] = ...
%             PerformBreathingRetrogatingAfterSliceReconForOnePhase(reconedSliceLocations, Nfe, Npe, NPhs, binAveAll, sliceThickness, ... 
%                         triggerTime, sliceLocation, cardiacPhases, binOverlapRatio, temporalRes, dataAcquired, option.resporitoryAcceptanceWindow, 0, ...
%                         strategy, inverse, initial, numOfPre, [100 64 1], 12, neighbor, stepDiv, moreIterInv, algo, volumePreserving, ...
%                         bestHeartBeatsImagePositionPatient, bestHeartBeatsImageOrientationPatient, mocoDir, selectedPhase);
                    
%% scatter interpolation
dstPixelSpacing = pixelSpacingRecon;
numOfSubdivision = 8;

[volumeBFFD, headerBFFD] = PerformScatterInterpolationFFD(binAveAll_BreathingGating_MocoAve, pixelSpacingRecon, ...
    binAveAll_BreathingGating_ImagePositionPatient, binAveAll_BreathingGating_ImageOrientationPatient, ... 
    dstPixelSpacing, numOfSubdivision, 'D:/temp/4DRecon/');

volumeBFFD = single(volumeBFFD);
save(fullfile(mocoDir, 'MoCOAveBFFD.mat'), 'volumeBFFD', 'headerBFFD');
% save4DVolume(volumeBFFD, headerBFFD, mocoDir, 'MocoAveBFFD', sizeRatio, centre, width, delay);

Matlab_SaveAnalyze(single(volumeBFFD), headerBFFD, fullfile(mocoDir, 'MocoAveBFFD.hdr'));

% numOfSubdivision = 8;
% [volumeLinearBFFD, headerLinearBFFD] = PerformScatterInterpolationFFD(binAveAll_BreathingGating_Linear, pixelSpacingRecon, ...
%     binAveAll_BreathingGating_ImagePositionPatient, binAveAll_BreathingGating_ImageOrientationPatient, ... 
%     dstPixelSpacing, numOfSubdivision, 'D:/temp/4DRecon/');
% 
% volumeLinearBFFD = single(volumeLinearBFFD);
% save(fullfile(mocoDir, 'LinearAveBFFD.mat'), 'volumeLinearBFFD', 'headerLinearBFFD');
% save4DVolume(volumeLinearBFFD, headerLinearBFFD, mocoDir, 'LinearBFFD', sizeRatio, centre, width, delay);

%% save the results
% save4DVolume(binAveAll_BreathingGating_MocoAve, header_binAveAll_BreathingGating, mocoDir, 'MocoAve', sizeRatio, centre, width, delay);
% save4DVolume(binAveAll_BreathingGating_Linear, header_binAveAll_BreathingGating, mocoDir, 'Linear', sizeRatio, centre, width, delay);
                    
% [startSLC, minInd] = min(reconedSliceLocations);
% [endSLC, maxInd] = max(reconedSliceLocations);
% 
% slices = [startSLC:(sign(endSLC-startSLC)*sliceThickness):endSLC];
% NSlc = numel(slices);
% 
% binAveAll_BreathingGating_Linear = zeros(Nfe, Npe, NPhs, NSlc);
% binAveAll_BreathingGating_MocoAve = zeros(Nfe, Npe, NPhs, NSlc);
% for s=1:NSlc
%     
%     if ( slices(s) <= startSLC )
%         binAveAll_BreathingGating(:,:,:,s) = binAveAll(minInd);
%         continue;
%     end
%     
%     if ( slices(s) >= endSLC )
%         binAveAll_BreathingGating(:,:,:,s) = binAveAll(maxInd);
%     end
%     
%     for k=1:numel(reconedSliceLocations)-1
%         if ( slices(s)>=reconedSliceLocations(k) & slices(s)<reconedSliceLocations(k+1) )
%             break;
%         end
%     end
%     
%     alpha = (slices(s)-reconedSliceLocations(k))/(reconedSliceLocations(k+1)-reconedSliceLocations(k));
%     binAveAll_BreathingGating_Linear(:,:,:,s) = (1-alpha)*binAveAll(:,:,:,k) + alpha*binAveAll(:,:,:,k+1);
%     
%     keyFrames = binAveAll_BreathingGating_Linear(:,:,:,s);
%     
%     % perform another round of retro-moco and averaging
%     currSliceRange = [slices(s)-sliceThickness slices(s)+sliceThickness];
%     for b=1:NPhs
%                
%         indPhs = find( ( triggerTime>=(cardiacPhases(b)-binOverlapRatio*temporalRes) ) ... 
%             & ( triggerTime<(cardiacPhases(b)+binOverlapRatio*temporalRes) ) ...
%             & (sliceLocation>=currSliceRange(1)) & (sliceLocation<currSliceRange(2)) );
%         
%         numOfIm = numel(indPhs);       
%         
%         disp(['b ' num2str(b) ' s ' num2str(s) ' numofIm ' num2str(numOfIm)]);
%         
%         if ( numOfIm > 1 )
% 
%             binData = zeros([size(data) numOfIm+1]);    
%             for ii=1:numOfIm
%                 dataX = dataAcquired(:,:,indPhs(ii));
%                 binData(:,:,ii) = dataX;
%             end            
%             
%             % get the keyframe image
%             currKeyFrame = keyFrames(:,:,b);
%             binData(:,:,end) = currKeyFrame;
%             
%             % perform the moco
%             header.sizeZ = size(binData, 3);
% 
%             [moco, dx, dy, dxInv, dyInv] = Matlab_PerformTemporalMotionCorrection(binData, header, numOfIm, strategy, inverse, initial, numOfPre, ...
%                 iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving);
% 
%             % perform the average
%             binAve = mean(moco(:,:,1:end-1), 3);
%             binAveAll_BreathingGating_MocoAve(:,:,b,s) = binAve;
%         end
%     end
% end
% 
% % switch the bin and slice
% binAveAll_BreathingGating_Linear = permute(binAveAll_BreathingGating_Linear, [1 2 4 3]);
% binAveAll_BreathingGating_MocoAve = permute(binAveAll_BreathingGating_MocoAve, [1 2 4 3]);
% 
% % compute the image position and orientation for the recon volume
% firstSlcImagePosition = bestHeartBeatsImagePositionPatient(minInd, :);
% firstSlcImageOrientation = bestHeartBeatsImageOrientationPatient(minInd, :);
% lastSlcImagePosition = bestHeartBeatsImagePositionPatient(maxInd, :);
% 
% header_binAveAll_BreathingGating = header;
% header_binAveAll_BreathingGating.sizeX = size(binAveAll_BreathingGating_MocoAve, 2);
% header_binAveAll_BreathingGating.sizeY = size(binAveAll_BreathingGating_MocoAve, 1);
% header_binAveAll_BreathingGating.sizeZ = size(binAveAll_BreathingGating_MocoAve, 3);
% header_binAveAll_BreathingGating.sizeT = size(binAveAll_BreathingGating_MocoAve, 4);
% header_binAveAll_BreathingGating.spacingZ = sliceThickness;
% header_binAveAll_BreathingGating.spacingT = temporalRes;
% header_binAveAll_BreathingGating.positionPatient = firstSlcImagePosition;
% header_binAveAll_BreathingGating.orientationPatient(:,1) = firstSlcImageOrientation(1:3); % row vector, x
% header_binAveAll_BreathingGating.orientationPatient(:,2) = firstSlcImageOrientation(4:6); % col vector, y
% % norm, z
% header_binAveAll_BreathingGating.orientationPatient(:,3) = (lastSlcImagePosition - firstSlcImagePosition) / norm(lastSlcImagePosition - firstSlcImagePosition); % row vector
% 
% filename = fullfile(mocoDir, ['binAveAll_BreathingGating_MocoAve']);
% Matlab_SaveAnalyze(single(binAveAll_BreathingGating_MocoAve), header_binAveAll_BreathingGating, [filename '.hdr']);
% 
% filename = fullfile(mocoDir, ['binAveAll_BreathingGating_Linear']);
% Matlab_SaveAnalyze(single(binAveAll_BreathingGating_Linear), header_binAveAll_BreathingGating, [filename '.hdr']);
% 
% Im4D = binAveAll_BreathingGating_MocoAve;
% save Im4D_BreathingGating_Linear Im4D header temporalRes binAveAll_BreathingGating_MocoAve binAveAll_BreathingGating_Linear

%% save results
% save binAveAll binAveAll
% save binAveAll_BreathingGating binAveAll_BreathingGating
% 
% binAveAllOri = binAveAll;
% binAveAll = binAveAll_BreathingGating_MocoAve;
% 
% for s=1:NSlc
%     binSlc = binAveAll(:,:,:,s);
%     header.sizeZ = size(binSlc, 3);
%     header.spacingZ = sliceThickness;
%     filename = fullfile(mocoDir, ['bin_SliceLocation' num2str(s) '_' suffix '_ave']);
%     Matlab_SaveAnalyze(single(binSlc), header, [filename '.hdr']);
%     hdr2gifWithWindowSetting([filename '.hdr'], [filename '.gif'], sizeRatio, centre, width, delay);
% end
% 
% for b=1:NPhs   
%     binAveAllTrigger = binAveAll(:,:,b,:);
%     binAveAllTrigger = squeeze(binAveAllTrigger);
%     
%     header.sizeZ = size(binAveAllTrigger, 3);
%     header.spacingZ = sliceThickness;
%     
%     filename = fullfile(mocoDir, ['bin_TriggerTime' num2str(b) '_' suffix '_ave']);
%     Matlab_SaveAnalyze(single(binAveAllTrigger), header, [filename '.hdr']);
%     hdr2gifWithWindowSetting([filename '.hdr'], [filename '.gif'], sizeRatio, centre, width, delay);
% end

%% a resporitory correction moco
% strategy = 'Consecutive';
% inverse = 1;
% initial = 0;
% numOfPre = 0;
% iters = [64 32 4];
% sigma = header.sizeY/4+header.sizeX/4;
% neighbor = 2.0;
% stepDiv = 3.0;
% moreIterInv = 1;
% algo = 'GLCC';
% volumePreserving = 1;
% 
% binAveAll_repoCorr = binAveAll;
% 
% % find the best slice as the key frames
% ssd = zeros(NSlc, NSlc);
% for s=1:NSlc    
%     binAveAllSlc = binAveAll(:,:,:,s);    
%     for s2=1:NSlc
%         if ( s ~= s2 )
%             binAveAllSlc2 = binAveAll(:,:,:,s2);
%             diffSlcIm = binAveAllSlc-binAveAllSlc2;
%             ssd(s, s2) = norm(diffSlcIm(:));
%         end
%     end        
% end
% 
% ssd = sort(ssd, 1);
% mssd = median(ssd);
% [minssd, bestSLC] = min(mssd);
% bestSLC = bestSLC - 1;
%     
% for b=1:NPhs
%     
%     disp(['Repo correction for trigger time ' num2str(b)]);
%     
%     binAveAllTrigger = binAveAll(:,:,b,:);
%     binAveAllTrigger = squeeze(binAveAllTrigger);
%     
% %     numOfIm = size(binAveAllTrigger, 3);
% %     ssd = zeros(numOfIm, numOfIm);
% %     for ff=1:numOfIm
% %         for ff2=1:numOfIm
% %             diffIm = binAveAllTrigger(:,:,ff)-binAveAllTrigger(:,:,ff2);
% %             ssd(ff, ff2) = norm(diffIm(:));
% %         end
% %     end
% % 
% %     ssd = sort(ssd, 1);
% %     mssd = median(ssd);
% %     [minssd, keyFrame] = min(mssd);
% %     keyFrame = keyFrame - 1;
% %     keyFrame
% 
%     header.sizeZ = size(binAveAllTrigger, 3);
%     
%     [moco, dx, dy, dxInv, dyInv] = Matlab_PerformTemporalMotionCorrection(double(binAveAllTrigger), header, bestSLC, strategy, inverse, initial, numOfPre, ...
%         iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving);
% 
%     binAveAll_repoCorr(:,:,b,:) = moco;
% end
% 
% %% save the results
% for s=1:NSlc   
%     binAveAllSlc = binAveAll_repoCorr(:,:,:,s);
%     
%     header.sizeZ = size(binAveAllSlc, 3);
%     header.spacingZ = sliceThickness;
%     
%     filename = fullfile(mocoDir, ['bin_Slice' num2str(s) '_' suffix '_ave_breathingCorrMoCo']);
%     Matlab_SaveAnalyze(single(binAveAllSlc), header, [filename '.hdr']);
%     hdr2gifWithWindowSetting([filename '.hdr'], [filename '.gif'], sizeRatio, centre, width, delay);
% end
% 
% for b=1:NPhs   
%     binAveAllTrigger = binAveAll_repoCorr(:,:,b,:);
%     binAveAllTrigger = squeeze(binAveAllTrigger);
%     
%     header.sizeZ = size(binAveAllTrigger, 3);
%     header.spacingZ = sliceThickness;
%     
%     filename = fullfile(mocoDir, ['bin_TriggerTime' num2str(b) '_' suffix '_ave_breathingCorrMoCo']);
%     Matlab_SaveAnalyze(single(binAveAllTrigger), header, [filename '.hdr']);
%     hdr2gifWithWindowSetting([filename '.hdr'], [filename '.gif'], sizeRatio, centre, width, delay);
% end
% 
% % binAveAll_repoCorr2 = binAveAll_repoCorr;
% % for s=1:NSlc
% %     
% %     binAveAllTrigger = binAveAll_repoCorr(:,:,:,s);
% %     binAveAllTrigger = squeeze(binAveAllTrigger);
% %     
% %     numOfIm = size(binAveAllTrigger, 3);
% %     ssd = zeros(numOfIm, numOfIm);
% %     for ff=1:numOfIm
% %         for ff2=1:numOfIm
% %             diffIm = binAveAllTrigger(:,:,ff)-binAveAllTrigger(:,:,ff2);
% %             ssd(ff, ff2) = norm(diffIm(:));
% %         end
% %     end
% % 
% %     ssd = sort(ssd, 1);
% %     mssd = median(ssd);
% %     [minssd, keyFrame] = min(mssd);
% %     keyFrame = keyFrame - 1;
% %     keyFrame
% % 
% %     header.sizeZ = size(binAveAllTrigger, 3);
% %     
% %     [moco, dx, dy, dxInv, dyInv] = Matlab_PerformTemporalMotionCorrection(double(binAveAllTrigger), header, keyFrame, strategy, inverse, initial, numOfPre, ...
% %         iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving);
% %     filename = fullfile(mocoDir, ['bin_SliceLocation' num2str(s) '_' suffix '_ave_breathingCorrMoCo']);
% %     Matlab_SaveAnalyze(single(moco), header, [filename '.hdr']);
% %     hdr2gifWithWindowSetting([filename '.hdr'], [filename '.gif'], sizeRatio, centre, width, delay);
% % 
% %     binAveAll_repoCorr2(:,:,:,s) = moco;
% % end
% % 
% % for b=1:NPhs
% %     binPhs = binAveAll_repoCorr2(:,:,b,:);
% %     header.sizeZ = size(binPhs, 3);
% %     filename = fullfile(mocoDir, ['bin_SliceLocation' num2str(s) '_' suffix '_ave_breathingCorrMoCo']);
% %     Matlab_SaveAnalyze(single(binSlc), header, [filename '.hdr']);
% %     hdr2gifWithWindowSetting([filename '.hdr'], [filename '.gif'], sizeRatio, centre, width, delay);
% % end
% 
% Im4D = binAveAll_repoCorr;
% 
% save binAveAll_repoCorr binAveAll_repoCorr

