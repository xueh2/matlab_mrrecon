function timeUsed = run_gt_recon_case(dataName, configName, UTCases, deleteh5)
% run the gt recon for a case

UTDir = getenv('GTPLUS_UT_DIR')

num = size(UTCases, 1);

for ii=1:num    
    ind = strfind(dataName, UTCases{ii, 2});
    if ( ~isempty(ind) & strcmp(configName, UTCases{ii, 4}) )
        testInfo = UTCases(ii, :);
        folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2})
        dataName = fullfile(folderDir, [UTCases{ii, 2} '.dat']);
        h5Name = fullfile(folderDir, [UTCases{ii, 2} '.h5']);
        resDir = fullfile(folderDir, UTCases{ii, 5});
        refDir = fullfile(folderDir, UTCases{ii, 6});
        xslFile = UTCases{ii, 7};

        isVD = strcmp(UTCases{ii, 3}, 'VD');
        isVD11 = strcmp(UTCases{ii, 3}, 'VD11');
        isAdjScan = strcmp(UTCases{ii, 3}, 'Adj');
        
        delete(fullfile(folderDir, 'res*.h5'));
        delete(fullfile(folderDir, 'out*.h5'));
        delete(fullfile(folderDir, '*.xml'));

        delete(fullfile(resDir, '*.nii'));
        delete(fullfile(resDir, 'gadgetron_*.hdr'));
        delete(fullfile(resDir, 'gadgetron_*.img'));
        delete(fullfile(resDir, '*.attrib'));

        if(deleteh5)
            delete(fullfile(resDir, '*.xml'));
            delete(fullfile(resDir, '*.h5'));
        end
                
        timeUsed = run_gt_recon(folderDir, dataName, h5Name, deleteh5, isVD, isVD11, isAdjScan, configName, resDir, xslFile, 1, 0, [], 0);
    end    
end
