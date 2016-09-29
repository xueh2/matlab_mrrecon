function run_gt_recon_VD13_All(folderName, configName, res, deleteh5)
% run the gt recon for a case

if nargin < 3
    res = 'res';
end

if nargin < 4
    deleteh5 = 0;
end

UTDir = getenv('GTPLUS_UT_DIR')

[names, num] = findFILE(folderName, '*.dat');

for ii=1:num
    dataName = names{ii};
    
    [path, name, ext] = fileparts(dataName);

    folderDir = fullfile(path, name);
    mkdir(folderDir);

    h5Name = fullfile(folderDir, [name '.h5']);
    resDir = fullfile(folderDir, res);

    xslFile = 'IsmrmrdParameterMap_Siemens.xsl';

    isVD = 1;
    isVD11 = 0;

    delete(fullfile(folderDir, 'res*.h5'));
    delete(fullfile(folderDir, 'out*.h5'));
    delete(fullfile(folderDir, '*.xml'));

    delete(fullfile(resDir, '*.nii'));
    delete(fullfile(resDir, '*.hdr'));
    delete(fullfile(resDir, '*.img'));
    delete(fullfile(resDir, '*.xml'));
    delete(fullfile(resDir, '*.h5'));

    run_gt_recon(folderDir, dataName, h5Name, deleteh5, isVD, isVD11, configName, resDir, xslFile);
end
