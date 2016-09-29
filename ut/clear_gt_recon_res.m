function clear_gt_recon_res(UTCases)
% clear the gt recon results

UTDir = getenv('GTPLUS_UT_DIR')

num = size(UTCases, 1);

for ii=1:num    
    testInfo = UTCases(ii, :);
    folderDir = fullfile(UTDir, UTCases{ii, 1}, UTCases{ii, 2})
    dataName = fullfile(folderDir, [UTCases{ii, 2} '.dat']);
    h5Name = fullfile(folderDir, [UTCases{ii, 2} '.h5']);
    resDir = fullfile(folderDir, UTCases{ii, 5});
    refDir = fullfile(folderDir, UTCases{ii, 6});

    isVD = strcmp(UTCases{ii, 3}, 'VD');

    delete(h5Name);
    
    delete(fullfile(folderDir, 'res*.h5'));
    delete(fullfile(folderDir, 'out*.h5'));
    delete(fullfile(folderDir, '*.xml'));

    delete(fullfile(resDir, '*.nii'));
    delete(fullfile(resDir, '*.hdr'));
    delete(fullfile(resDir, '*.img'));
    delete(fullfile(resDir, '*.xml'));
    delete(fullfile(resDir, '*.h5'));
end
