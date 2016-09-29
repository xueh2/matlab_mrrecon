
function removeUTValidation(UTCases, res)
%% remove a UT validation folder

GTHome = getenv('GADGETRON_HOME')
addpath(fullfile(GTHome, 'bin'));
addpath(fullfile(GTHome, 'lib'));

UTDir = getenv('GTPLUS_UT_DIR')
addpath(UTDir)

%% perform the unit test

numOfCases = size(UTCases, 1);

% set the ref
startCase = 1;
endCase = numOfCases;

% dos('del c:\temp\gadgetron\*');

for ii=startCase:endCase

    folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2})
    resDir = fullfile(folderDir, res);
    
    try
        rmdir(resDir, 's');
    catch
    end
end
