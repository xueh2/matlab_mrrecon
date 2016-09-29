
clear all
close all

%% compare the unit test results with ground truth

GTHome = getenv('GADGETRON_HOME')
addpath(fullfile(GTHome, 'bin'));
addpath(fullfile(GTHome, 'lib'));

UTDir = getenv('GTPLUS_UT_DIR')
addpath(UTDir)

c = clock;
logfile = fullfile(UTDir, 'log', ['GtPlus_UT_' date '_' num2str(c(4)) '_' num2str(c(5)) '_' num2str(c(6)) '.txt']);
logfile
if( isFileExist(logfile) )
    command = ['del /F /Q ' logfile]; 
    delete(logfile);
end

setenv('OutputFormat', 'hdr');

% setenv('GT_HOST', 'localhost');
% setenv('GT_HOST', '137.187.134.9');
% setenv('GT_HOST', 'barbados.nhlbi.nih.gov');
% setenv('GT_HOST', 'nauru.nhlbi.nih.gov');
% setenv('GT_HOST', '128.231.26.242');

% setenv('GT_PORT', '9002');
% setenv('GT_PORT', '9004');

% setenv('GT_PORT', '9016');
% setenv('GT_PORT', '9222');

%% set the ut cases
% UTCases = set_up_UT_cases
UTCases = set_up_UT_cases_ScannerUpdate

%% perform the unit test

numOfCases = size(UTCases, 1);

% set the ref
UTCollectRef = 0;

for ii=1:numOfCases
    
    testInfo = UTCases(ii, :);
    folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2})
    dataName = fullfile(folderDir, [UTCases{ii, 2} '.dat']);
    h5Name = fullfile(folderDir, [UTCases{ii, 2} '.h5']);
    resDir = fullfile(folderDir, UTCases{ii, 5});
    refDir = fullfile(folderDir, UTCases{ii, 6});

    isVD = strcmp(UTCases{ii, 3}, 'VD');
    configName = UTCases{ii, 4};
        
    run_gt_recon_case(UTCases{ii, 2}, UTCases{ii, 4}, UTCases); % winopen(folderDir)
     
    if ( UTCollectRef )
        collect_UT_reference(UTCollectRef, resDir, refDir);
    else
        compare_with_reference(testInfo, resDir, refDir, logfile);
    end
end
winopen(logfile);

error('stop here ...');

data = readGTPlusExportImages('.', 'GT_ImageRetro'); data=squeeze(data); figure; imagescn(data, [], [], [], 3);
data = readGTPlusExportImages('.', 'GT_ImageRetro'); data=squeeze(data); figure; imagescn(data, [], [], [], 4);

RR = 1000;
TR = 97;
alternating = 0;
pmuTime = KSpaceBinningSimu(20000, RR, TR, 110, 5, alternating);
ind = find(pmuTime(:)>0);
% plot(pmuTime(ind(:)), '.')

dstPHS = 30;
fillingMatrix = zeros(size(pmuTime, 1), dstPHS);
for phs=1:dstPHS
    IND = find(pmuTime(:)<phs*RR/dstPHS & pmuTime(:)>(phs-1)*RR/dstPHS);
    [I, J] = ind2sub(size(pmuTime), IND);
    
    N = numel(I);
    for n=1:N
        fillingMatrix(I(n), phs) = fillingMatrix(I(n), phs) + 1;
    end
end

imtool(fillingMatrix, 'InitialMagnification', 400)

