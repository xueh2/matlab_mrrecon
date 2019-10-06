
function performUTMD5SUM(UTCases, caseNum, caseNumEnd)
%% compute md5 sum for ut results
% performUTMD5SUM(UTCases, caseNum, caseNumEnd)

GTHome = getenv('GADGETRON_HOME')
UTDir = getenv('GTPLUS_UT_DIR')
OutputFormat = getenv('OutputFormat');

if nargin < 2
    caseNum = 1;
end

if nargin < 3
    caseNumEnd = size(UTCases, 1);
end

%% perform the unit test

numOfCases = size(UTCases, 1);

% set the ref
startCase = 1;
endCase = numOfCases;

if ( caseNum > 0 )
    startCase = caseNum;
    endCase = caseNumEnd;
end

date_suffix = datestr(date, 'yyyymmdd');

for ii=startCase:endCase

    testInfo = UTCases(ii, :);
    folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2})
    
    if(~isFileExist(folderDir))
        folderDir = fullfile(UTDir, UTCases{ii, 1});
    end
    
    dataName = fullfile(folderDir, [UTCases{ii, 2} '.dat']);
    resDir = fullfile(folderDir, [UTCases{ii, 5}]);

    [names, num] = findFILE(resDir, ['*' date_suffix '*.h5']);
    
    ind = find(folderDir=='\');
    folderDir(ind) = '/';
    
    for k=1:num
        disp(['=====================================================']);
        disp(folderDir)
        
        ind = find(names{k}=='\');
        fname = names{k};
        fname(ind) = '/';
        dos(['fciv.exe -md5 ' fname], '-echo');
    end
end
