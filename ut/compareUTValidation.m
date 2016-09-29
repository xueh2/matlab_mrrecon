
function [a, b] = compareUTValidation(UTCases, res1, res2, seriesNo, timeDim, windowing, remove_scaling, startCase, endCase)
% compare the unit test results with ground truth
% compareUTValidation(set_up_UT_cases_ImageReviewing, 'res', 'res_old', seriesNo, timeDim, windowing, remove_scaling, startCase, endCase)

%% run
GTHome = getenv('GADGETRON_HOME')
addpath(fullfile(GTHome, 'bin'));
addpath(fullfile(GTHome, 'lib'));

UTDir = getenv('GTPLUS_UT_DIR')
addpath(UTDir)

%% perform the unit test

numOfCases = size(UTCases, 1);

if(nargin<7)
    remove_scaling = 0;
end

if(nargin<8)
    startCase = 1;
end

if(nargin<9)
    endCase = numOfCases;
end

% dos('del c:\temp\gadgetron\*');

timeUsed = cell(numOfCases, 4);

for ii=startCase:endCase

    testInfo = UTCases(ii, :);
    folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2})
    dataName = fullfile(folderDir, [UTCases{ii, 2} '.dat']);
    h5Name = fullfile(folderDir, [UTCases{ii, 2} '.h5']);
    resDir1 = fullfile(folderDir, res1);
    resDir2 = fullfile(folderDir, res2);
    
    a = [];
    b = [];
    if(isdir(resDir1))
        a = readGTPlusExportImageSeries_Squeeze(resDir1, seriesNo);
    end
    
    if(isdir(resDir2))
        b = readGTPlusExportImageSeries_Squeeze(resDir2, seriesNo);
    end
    
    if(isempty(a))
        figure; imagescn(b, [windowing], [], [], timeDim);
    elseif (isempty(b))
        figure; imagescn(a, [windowing], [], [], timeDim);
    else
        if(remove_scaling)
            sa = norm(a(:));
            sb = norm(b(:));
            b = b * sa/sb;
        end
        figure; imagescn(cat(numel(size(a))+1, a, b, (windowing(2)+windowing(1))/2 + (a-b)*300 ), [windowing], [], 1600/size(a, 1), timeDim);
    end
    pause
    try
        [I, xlim, ylim, clim, position]=getimagedata;
        windowing = [clim(1, 1) clim(1, 2)];
    catch
    end
    closeall
end
