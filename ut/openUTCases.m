function openUTCases(UTCases, casename)
% open the UT cases folder
% openUTCases(set_up_UT_cases_ScannerUpdate, casename)

UTDir = getenv('GTPLUS_UT_DIR')
addpath(UTDir)

numOfCases = size(UTCases, 1);

for ii=1:numOfCases
    
    if ( strcmp(casename, UTCases{ii, 2}) == 1 )       
        winopen(fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2}, UTCases{ii, 5}));       
        break;
    end
    
end