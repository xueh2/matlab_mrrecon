function run_gt_recon_case_ICE(dataName, UTCases)
% run the gt recon for a case

UTDir = getenv('GTPLUS_UT_DIR')
if ( numel(UTDir) == 0 )
    UTDir = 'E:/ut/';
end

num = size(UTCases, 1);

isVD11 = 0;

for ii=1:num    
    ind = strfind(dataName, UTCases{ii, 2});
    if ( ~isempty(ind) )
        testInfo = UTCases(ii, :);
        folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2})
        dataName = fullfile(folderDir, [UTCases{ii, 2} '.dat']);
        resDir = fullfile(folderDir, UTCases{ii, 3});
        mkdir(resDir);
        refDir = fullfile(folderDir, UTCases{ii, 4});

        delete(fullfile(resDir, '*.hdr'));
        delete(fullfile(resDir, '*.img'));
        delete(fullfile(resDir, '*.xml'));
        delete(fullfile(resDir, '*.dcm'));
        delete(fullfile(resDir, '*.Icehead'));
        delete(fullfile(resDir, '*.ima'));

        cd('Z:\n4\pkg\MrServers\MrVista\Ice\temp')
        command = ['icesimu ' dataName];
        dos(command, '-echo');
        
        collect_UT_reference(1, folderDir, resDir);        
        
        delete(fullfile(folderDir, '*.hdr'));
        delete(fullfile(folderDir, '*.img'));
        delete(fullfile(folderDir, '*.xml'));
        delete(fullfile(folderDir, '*.dcm'));
        delete(fullfile(folderDir, '*.Icehead'));
        delete(fullfile(folderDir, '*.ima'));

    end    
end
