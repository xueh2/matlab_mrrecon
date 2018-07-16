
function timeUsed = PerformGadgetronRecon_SavedMeasDat_OneData(dataDir, data_name, gt_host, resDir, configName, isPerf, styleSheet)
% timeUsed = PerformGadgetronRecon_SavedMeasDat_OneData(dataDir, data_name, gt_host, resDir, configName, isPerf, styleSheet)
% timeUsed = PerformGadgetronRecon_SavedIsmrmrd_OneData('I:\KAROLINSKA', 'DB_LGE_MOCO_AVE_OnTheFly_41672_7432041_7432046_590_20160330-113443', 'localhost', 'I:\ReconResults\KAROLINSKA')

GT_PORT = gtPortLookup(gt_host);
[key, user] = sshKeyLookup(gt_host);

if(nargin<4)
    resDir = dataDir;
end

if(nargin<6)
    isPerf = 0;
end

if(nargin<7)
    styleSheet = 'IsmrmrdParameterMap_Siemens_Perfusion.xsl';
end

setenv('GT_HOST', gt_host)
setenv('GT_PORT', GT_PORT)

GTHome = getenv('GADGETRON_HOME')
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

styleSheetDefault = 'IsmrmrdParameterMap_Siemens.xsl';
styleSheetPerfusionUsed = 'IsmrmrdParameterMap_Siemens_Perfusion.xsl';
if ( nargin >= 6 )
    styleSheetDefault = [ styleSheet];
end

xmlUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xml';
   
% ------------------------------------------------------------
  
if(isempty(data_name))
    [names, num] = findFILE(dataDir, '*.dat');
else
    [path, dat_name, ext] = fileparts(data_name);
    dataName = fullfile(dataDir, [data_name '.dat']);
    names{1} = dataName;
end

for ii=1:numel(names)
    [path, dat_name, ext] = fileparts(names{ii});
    dataName = names{ii};
    
    dstDir = fullfile(resDir, dat_name);
    mkdir(dstDir);
    cd(dstDir)

    delete(fullfile(dstDir, 'res*.h5'));
    delete(fullfile(dstDir, 'out*.h5'));
    delete(fullfile(dstDir, '*.xml'));

    delete(fullfile(dstDir, '*.nii'));
    delete(fullfile(dstDir, 'gadgetron_*.hdr'));
    delete(fullfile(dstDir, 'gadgetron_*.img'));
    delete(fullfile(dstDir, 'Generic*.hdr'));
    delete(fullfile(dstDir, 'Generic*.img'));
    delete(fullfile(dstDir, 'GTPrep*.hdr'));
    delete(fullfile(dstDir, 'GTPrep*.img'));
    delete(fullfile(dstDir, 'GT*.hdr'));
    delete(fullfile(dstDir, 'GT*.img'));
    delete(fullfile(dstDir, '*.attrib'));

    delete(fullfile(dstDir, '*.xml'));          

    finfo = dir(dataName);

    %% run the scan
    configNameUsed = fullfile(GTConfigFolder, configName);
    cd(dstDir)

    h5Name = fullfile(dstDir, [dat_name '.h5']);
    deleteh5 = 1;
    isVD = 1;
    isVD11 = 0;
    isAdjScan = 0;
    startRemoteGT = 0;
    h5Only = 0;
    remoteXml = [];
    compressionBit = [];
    
    if(isPerf)
        if(strcmp(getenv('GT_HOST'), 'localhost')==1)
            debugFolder = 'D:\gtuser\mrprogs\install\DebugOutput';
            try
                rmdir(debugFolder, 's');
            catch
            end

            try
                mkdir(debugFolder);
            catch
            end
        end
    end
    
    timeUsed = run_gt_recon(dstDir, dataName, h5Name, deleteh5, isVD, isVD11, isAdjScan, configName, dstDir, styleSheetDefault, startRemoteGT, h5Only, remoteXml, compressionBit);

    if(isPerf)
        if(strcmp(getenv('GT_HOST'), 'localhost')==1)
            movefile(debugFolder, dstDir, 'f');
        else
            [key, user] = sshKeyLookup(getenv('GT_HOST'));
            debug_folder = ['/home/' user '/Debug/DebugOutput']
            CopyGadgetronDebugOutputOnRemote(getenv('GT_HOST'), debug_folder, dstDir, 1)
        end
    end

    %% copy the dicom
    PerformGadgetronRecon_SavedIsmrmrd_CopyDicom(resDir, data_name, gt_host);

    %% open the results
    winopen(dstDir)
end
