
function [timeUsed, remoteFolder] = PerformGadgetronRecon_SavedIsmrmrd_CopyDicom(resDir, data_name, gt_host, remote_dicom_dir)
% timeUsed = PerformGadgetronRecon_SavedIsmrmrd_CopyDicom(resDir, data_name, gt_host, remote_dicom_dir)
% timeUsed = PerformGadgetronRecon_SavedIsmrmrd_CopyDicom('I:\ReconResults\KAROLINSKA', 'DB_LGE_MOCO_AVE_OnTheFly_41672_7432041_7432046_590_20160330-113443', 'barbados')
% timeUsed = PerformGadgetronRecon_SavedIsmrmrd_CopyDicom('I:\ReconResults\KAROLINSKA', 'DB_LGE_MOCO_AVE_OnTheFly_41672_7432041_7432046_590_20160330-113443', 'barbados', '/tmp/gadgetron_data')

[key, user] = sshKeyLookup(gt_host);

if(nargin<4)
    if(strcmp(gt_host, 'localhost')==1)
        remote_dicom_dir = 'd:/temp/gadgetron_data';
    else
        remote_dicom_dir = '/tmp/gadgetron_data';
    end
end

% ------------------------------------------------------------

[configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(data_name);

dstDir = fullfile(resDir, study_dates, [data_name '_dicom']);
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

%% copy the dicom

remoteFolder = data_name;
if(strcmp(gt_host, 'localhost')==1)
    tic; copyfile(fullfile(remote_dicom_dir, remoteFolder, '*.dcm'),  '.'); timeUsed = toc;    
else
    command = ['pscp -i ' key '.ppk ' user '@' gt_host ':' remote_dicom_dir '/' remoteFolder '/*.dcm .'];
    command
    tic; dos(command, '-echo'); timeUsed = toc;
end

dstDir

%% 
% winopen(dstDir)