
function timeUsed = performUTValidation(UTCases, UTCollectRef, deleteh5, GT_HOST, GT_PORT, caseNum, caseNumEnd, compressionBit, compRef, startRemoteGT, remoteXml, res_suffix, paraXml, debug_folder);
%% compare the unit test results with ground truth
% performUTValidation(UTCases, UTCollectRef, deleteh5, GT_HOST, GT_PORT, caseNum, caseNumEnd, compressionBit, compRef, startRemoteGT, remoteXml, res_suffix)
% performUTValidation(UTCases, UTCollectRef, GT_HOST, GT_PORT, caseNum)
% performUTValidation(set_up_UT_cases_ScannerUpdate, 0, 0, 'localhost', '9002')
% performUTValidation(set_up_UT_cases_ScannerUpdate, 0, 0, 'barbados.nhlbi.nih.gov', '9006')
% performUTValidation(set_up_UT_cases_Cloud, 0, 0, 'nauru.nhlbi.nih.gov', '9002')
% performUTValidation(set_up_UT_cases_GT_UT, 0, 0, 'localhost', '9002')
% performUTValidation(set_up_UT_cases_Cine, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0)
% performUTValidation(set_up_UT_cases_Perfusion, 0, 0, '137.187.135.157', '9008', 63, 63, 0, 0, 0)
% performUTValidation(UTCases, UTCollectRef, deleteh5, GT_HOST, GT_PORT, caseNum, caseNumEnd, compressionBit, compRef, startRemoteGT, remoteXml, res_suffix)
% performUTValidation(UTCases, UTCollectRef, deleteh5, GT_HOST, GT_PORT, caseNum, caseNumEnd, compressionBit, compRef, startRemoteGT, remoteXml, res_suffix, paraXml)

GTHome = getenv('GADGETRON_HOME')
% addpath(fullfile(GTHome, 'bin'));
% addpath(fullfile(GTHome, 'lib'));

UTDir = getenv('GTPLUS_UT_DIR')
% addpath(UTDir)

c = clock;
logfile = fullfile(UTDir, 'log', ['GtPlus_UT_' date '_' num2str(c(4)) '_' num2str(c(5)) '_' num2str(c(6)) '.txt']);
logfile
if( isFileExist(logfile) )
    command = ['del /F /Q ' logfile]; 
    delete(logfile);
end

% setenv('OutputFormat', 'hdr');

if nargin < 1
    UTCases = set_up_UT_cases_ScannerUpdate;
    UTCollectRef = 0;
    deleteh5 = 0;
    GT_HOST='localhost';
    GT_PORT='9002';
end

if nargin < 2
    UTCollectRef = 0;
    deleteh5 = 0;
    GT_HOST='localhost';
    GT_PORT='9002';
end

if nargin < 6
    caseNum = -1;
end

if nargin < 7
    caseNumEnd = caseNum;
end

if nargin < 8
    compressionBit = 0;
end

if nargin < 9
    compRef = 1;
end

if nargin < 10
    startRemoteGT = 1;
end

if nargin < 11
    remoteXml = 0;
end

if nargin < 12
    res_suffix = [];
end

if nargin < 13
    paraXml = [];
end

if nargin < 14
    debug_folder = [];
end

if(strcmp(GT_HOST, 'gt1'))
    GT_HOST = '137.187.135.169';
end
if(strcmp(GT_HOST, 'gt2'))
    GT_HOST = '137.187.135.238';
end
if(strcmp(GT_HOST, 'gt3'))
    GT_HOST = '137.187.135.194';
end
if(strcmp(GT_HOST, 'beast'))
    GT_HOST = '137.187.135.157';
end

setenv('GT_HOST', GT_HOST);
setenv('GT_PORT', GT_PORT);


%% perform the unit test

numOfCases = size(UTCases, 1);

% set the ref
startCase = 1;
endCase = numOfCases;

if ( caseNum > 0 )
    startCase = caseNum;
    endCase = caseNumEnd;
end

% dos('del c:\temp\gadgetron\*');

timeUsed = cell(numOfCases, 4);

for ii=startCase:endCase

    testInfo = UTCases(ii, :);
    folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2})
    
    if(~isFileExist(folderDir))
        folderDir = fullfile(UTDir, UTCases{ii, 1})
    end
    
    if(~isFileExist(folderDir))
        folderDir = fullfile(UTCases{ii, 1}, UTCases{ii, 2})
    end
    
    dataName = fullfile(folderDir, [UTCases{ii, 2} '.dat']);
    h5Name = fullfile(folderDir, [UTCases{ii, 2} res_suffix '.h5']);
    resDir = fullfile(folderDir, [UTCases{ii, 5} res_suffix]);
    refDir = fullfile(folderDir, [UTCases{ii, 6} res_suffix]);

    isVD = strcmp(UTCases{ii, 3}, 'VD');
    isAdj = strcmp(UTCases{ii, 3}, 'Adj');
    configName = UTCases{ii, 4};

    h5Only = 0;
    if(~isFileExist(dataName))
        h5Only = 1;
    end
    
    UT_case_Dir = [];
    t = run_gt_recon_case(UTCases{ii, 2}, UTCases{ii, 4}, UTCases, deleteh5, startRemoteGT, res_suffix, h5Only, remoteXml, compressionBit, paraXml, UT_case_Dir, debug_folder); % winopen(folderDir)
    disp(['--> run time ' num2str(t) ' seconds ... ']);
    
    if(isAdj) 
        continue; 
    end
    
    timeUsed{ii, 1} = ii;
    timeUsed{ii, 2} = t;
    timeUsed{ii, 3} = UTCases{ii, 2};
    timeUsed{ii, 4} = folderDir;
    
    if ( UTCollectRef )
        try
            collect_UT_reference(UTCollectRef, resDir, refDir);
        catch
        end
    else
        try
            if(compRef)
                compare_with_reference(testInfo, resDir, refDir, logfile);
            end
        catch
        end
    end
    
    if ( endCase==startCase )
        if(isunix())
            disp(resDir);
            system(['nautilus ' resDir])
        else
            winopen(resDir);
        end
    end
end

if ( UTCollectRef == 0 & compRef==1 )
    try
        dos(['code ' logfile]);
    catch
    end
end
