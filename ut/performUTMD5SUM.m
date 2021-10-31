
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

    disp(['=====================================================================================']);
    
    testInfo = UTCases(ii, :);
    folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2});
    
    if(~isFileExist(folderDir))
        folderDir = fullfile(UTDir, UTCases{ii, 1});
    end
    
    dataName = fullfile(folderDir, [UTCases{ii, 2} '.dat']);
    resDir = fullfile(folderDir, [UTCases{ii, 5}]);

    [names, num] = findFILE(resDir, ['*' date_suffix '*.h5']);
    
    ind = find(folderDir=='\');
    folderDir(ind) = '/';
    
    for k=1:num
        disp(['-------------------']);
        
        ind = find(names{k}=='\');
        fname = names{k};
        fname(ind) = '/';
        % dos(['fciv.exe -md5 ' fname], '-echo');
        if(isunix())
            [status, outputs] = system(['md5sum ' fname]);
        else
            %[status, outputs] = system(['fciv.exe -md5 ' fname]);
            [status, outputs] = system(['D:\software\md5sums.exe -u ' fname]);
        end
        print_json_record(outputs, fname);
    end
    
    disp(['******************']);
    
    if(exist(dataName))
        if(isunix())
            [status, outputs] = system(['md5sum ' dataName]);
        else
            %[status, outputs] = system(['fciv.exe -md5 ' dataName]);
            [status, outputs] = system(['D:\software\md5sums.exe -u ' dataName]);
        end
        print_json_record(outputs, dataName);
    else
        if(isunix())
            [status, outputs] = system(['md5sum ' fullfile(folderDir, [UTCases{ii, 2} '.h5'])]);
        else
            %[status, outputs] = system(['fciv.exe -md5 ' fullfile(folderDir, [UTCases{ii, 2} '.h5'])]);
            [status, outputs] = system(['D:\software\md5sums.exe -u ' fullfile(folderDir, [UTCases{ii, 2} '.h5'])]);
        end
        print_json_record(outputs, fullfile(folderDir, [UTCases{ii, 2} '.h5']));
    end
    
    
end

end

function print_json_record(outputs, fname)
    ind = strfind(outputs, '//');
    if(isempty(ind))
        ind = strfind(outputs, '*');
        if(isempty(ind))
            ind = strfind(outputs, '/mnt');
            a = outputs(ind(end):end);
            ind2 = strfind(outputs, ' ');
            sha1 = outputs(1:ind2(1)-1);
        else
            a = outputs(1:ind(1)-1);
            ind2 = strfind(a, ' ');
            sha1 = a(1:ind2(1)-1);
        end
    else
        a = outputs(ind(end)+3:end);
        ind2 = strfind(a, ' ');
        sha1 = a(1:ind2(1)-1);
    end
    ind3 = strfind(fname, 'gtprep/');
    if(isempty(ind3))
        ind3 = strfind(fname, 'gtprep\');
    end
    fname_str = fname(ind3:end);
    ind = strfind(fname_str, '\');
    fname_str(ind) = '/';
    disp(['{']);
    disp(['    "file":"' fname_str '",']);
    disp(['    "md5":"' sha1 '"']);
    disp(['}']);
end
