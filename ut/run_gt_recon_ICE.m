function run_gt_recon_ICE(dataName, deleteh5)
% run the gt recon for a case

if  nargin == 1
    deleteh5 = 0;
end

UTDir = getenv('GTPLUS_UT_DIR')
if ( numel(UTDir) == 0 )
    UTDir = 'D:/gtuser/gt_windows_setup/ut/';
end

styleSheetUsed = '%GADGETRON_DIR%\install\gadgetron\schema/IsmrmrdParameterMap_Siemens.xsl';

[folderDir, name, ext] = fileparts(dataName);

resDir = fullfile(folderDir, 'ICE_res');
mkdir(resDir);

delete(fullfile(resDir, '*.hdr'));
delete(fullfile(resDir, '*.img'));
delete(fullfile(resDir, '*.xml'));
delete(fullfile(resDir, '*.dcm'));
delete(fullfile(resDir, '*.Icehead'));
delete(fullfile(resDir, '*.ima'));

h5Name = fullfile(folderDir, [name '.h5']);
if ( ~isFileExist(h5Name) || deleteh5 )
    delete(h5Name);
    command = ['siemens_to_HDF5 ' dataName ' ' h5Name]
    dos(command);
end

command = ['siemens_mriclient -f ' h5Name ' -m %GADGETRON_DIR%\install\gadgetron\schema/IsmrmrdParameterMap_Siemens.xml -x %GADGETRON_DIR%\install\gadgetron\schema/IsmrmrdParameterMap_Siemens.xsl -d 0 -c default_measurement_dependencies.xml -h %GT_HOST% -p %GT_PORT% ']
dos(command);
        
dos('icestop');
dos('icestart -r -dicom');

cd('Y:\n4\pkg\MrServers\MrVista\Ice\temp')
command = ['icesimu ' dataName];
dos(command, '-echo');

collect_UT_reference(1, folderDir, resDir);        

delete(fullfile(folderDir, '*.hdr'));
delete(fullfile(folderDir, '*.img'));
delete(fullfile(folderDir, '*.xml'));
delete(fullfile(folderDir, '*.dcm'));
delete(fullfile(folderDir, '*.Icehead'));
delete(fullfile(folderDir, '*.ima'));

dos('icestop');

winopen(resDir);

