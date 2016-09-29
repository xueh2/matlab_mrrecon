
function cardiacPhases = DetermineJointCardiacPhase_20120429(dataDirs, temporalRes, plotFlag)
% function cardiacPhases = DetermineJointCardiacPhase_20120429(dataDirs)
%
% this function compute the joint mean destination cardiac phase for multiple acquisitions
%
% Inputs:
%    dataDirs : folders storing the 4D datasets in dicom format (*.ima or *.dcm files)
%   temporalRes : required temproal resolution
%
% Output:
%    cardiacPhases : [PHS 1] vector, storing the cardiac phases reconstructed
%
%     ***************************************
%     *  Hui Xue (hui-xue@siemens.com       *
%     *  2012-03                            *
%     ***************************************

%% parameters

N = numel(dataDirs);

startingPhaseAll = [];
endingPhaseAll = [];

for n=1:N
    
dataDir = dataDirs{n};
cd(dataDir)
currDir = dataDir;

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
    save data4D_TriggerTime triggerTime
else
    load(filename);
end


%% compute the mean R-R interval, mean starting phase and mean ending phase and other cardiac phases
disp('-----------------------------------------');
% all indexes are in the order of global frames
[meanRR, meanTriggerTimeDelta, startingPhase, startingInd, endingPhase, endingInd, indexes] = ComputeMeanRR(triggerTime);

save CardiacPhase startingPhase endingPhase

startingPhaseAll = [startingPhaseAll; startingPhase];
endingPhaseAll = [endingPhaseAll; endingPhase];

end

% compute the mean cardiac phase
startPhase = mean(startingPhaseAll);
endPhase = mean(endingPhaseAll);

cardiacPhases = [startPhase:temporalRes:endPhase];
NPhs = numel(cardiacPhases);
disp(['Number of reconstructed cardiac phases : ' num2str(numel(cardiacPhases))]);
