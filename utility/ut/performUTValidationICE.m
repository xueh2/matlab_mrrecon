
function performUTValidationICE(UTCases, UTCollectRef, caseNum, caseEnd)
% compare the unit test results with ground truth
% performUTValidationICE(set_up_UT_cases_ICE, UTCollectRef, caseNum)

GTHome = getenv('GADGETRON_HOME')
addpath(fullfile(GTHome, 'bin'));
addpath(fullfile(GTHome, 'lib'));

UTDir = getenv('GTPLUS_UT_DIR')
if ( numel(UTDir) == 0 )
    UTDir = 'E:/ut/';
end
addpath(UTDir)

c = clock;
logfile = fullfile(UTDir, 'log', ['GtPlus_UT_' date '_' num2str(c(4)) '_' num2str(c(5)) '_' num2str(c(6)) '.txt']);
logfile
if( isFileExist(logfile) )
    command = ['del /F /Q ' logfile]; 
    delete(logfile);
end

if nargin < 1
    UTCases = set_up_UT_cases_ICE;
    UTCollectRef = 0;
end

if nargin < 2
    UTCollectRef = 0;
end

if nargin < 3
    caseNum = -1;
end

if nargin < 4
    caseEnd = -1;
end

%% perform the unit test

numOfCases = size(UTCases, 1);

% set the ref
startCase = 1;
endCase = numOfCases;

if ( caseNum > 0 )
    startCase = caseNum;
end

if ( caseEnd > 0 )
    endCase = caseEnd;
end

dos('icestop', '-echo');
dos('icestart -r -dicom', '-echo');

for ii=startCase:endCase

    testInfo = UTCases(ii, :);
    folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2})
    dataName = fullfile(folderDir, [UTCases{ii, 2} '.dat']);
    resDir = fullfile(folderDir, UTCases{ii, 3});
    refDir = fullfile(folderDir, UTCases{ii, 4});

    run_gt_recon_case_ICE(UTCases{ii, 2}, UTCases); % winopen(folderDir)

    if ( UTCollectRef )
        collect_UT_reference(UTCollectRef, resDir, refDir);
    else
        try
            compare_with_reference_ICE(testInfo, resDir, refDir, logfile);
        catch
        end
    end
    
    if ( caseNum > 0 )
        winopen(resDir);
    end
end

dos('icestop', '-echo');

if ( UTCollectRef == 0 )
    winopen(logfile);
end
