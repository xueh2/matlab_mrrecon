
function md5 = copy_for_UT_ref(UTCases, dstDir)
%% copy results for ut ref
% copy_for_UT_ref(UTCases, dstDir)
% copy_for_UT_ref(set_up_UT_cases_GtPrep_UT, 'D:\data\gtprep_ut_20171116')

GTHome = getenv('GADGETRON_HOME')
UTDir = getenv('GTPLUS_UT_DIR')

numOfCases = size(UTCases, 1);

% set the ref
startCase = 1;
endCase = numOfCases;

md5 = [];

for ii=startCase:endCase

    testInfo = UTCases(ii, :);
    folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2})
    
    if(~isFileExist(folderDir))
        folderDir = fullfile(UTDir, UTCases{ii, 1})
    end
    
    dataName = fullfile(folderDir, [UTCases{ii, 2} '.dat']);
    h5Name = fullfile(folderDir, [UTCases{ii, 2} '.h5']);
    resDir = fullfile(folderDir, [UTCases{ii, 5} ]);
    refDir = fullfile(folderDir, [UTCases{ii, 6} ]);

    dstDir_case = fullfile(dstDir, UTCases{ii, 1}, UTCases{ii, 2});
    mkdir(dstDir_case);
    mkdir(fullfile(dstDir_case, 'ref'));
    
    try
        copyfile(dataName, fullfile(dstDir_case, [UTCases{ii, 2} '.dat']));
        has_dat = 1;
    catch
        copyfile(h5Name, fullfile(dstDir_case, [UTCases{ii, 2} '.h5']));
        has_dat = 0;
    end
    copyfile(fullfile(resDir, 'ref*.h5'), fullfile(dstDir_case, 'ref'));
    
    a1 = mMD5(fullfile(dstDir_case, [UTCases{ii, 2} '.dat']));
%     command = ['md5sum ' fullfile(dstDir_case, [UTCases{ii, 2} '.dat'])];
%     dos(command, '-echo');
    
    [names, num] = findFILE(fullfile(dstDir_case, 'ref'), '*.h5');
    b1 = mMD5(names{1});
    
    [path, name, ext] = fileparts(names{1});
    
    if(has_dat)
        m_str = [UTCases{ii, 1} '/' UTCases{ii, 2} '/' UTCases{ii, 2} '.dat: ' a1];
    else
        m_str = [UTCases{ii, 1} '/' UTCases{ii, 2} '/' UTCases{ii, 2} '.h5: ' a1];
    end
    ind = find(m_str=='\');
    m_str(ind) = '/';
    
    md5 = [md5; {m_str}];
    
    m_str = [UTCases{ii, 1} '/' UTCases{ii, 2} '/ref/' name '.h5: ' b1];
    ind = find(m_str=='\');
    m_str(ind) = '/';
    md5 = [md5; {m_str}];
%     command = ['md5sum ' names{1}];
%     dos(command, '-echo');

    
end
