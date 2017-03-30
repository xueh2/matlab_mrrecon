
function PerfCases = PerformGadgetronRecon_FindPerfusion_PatientID_PerfusionCase(dataDir, resDir, patientID, study_date)
% PerfCases = PerformGadgetronRecon_FindPerfusion_PatientID_PerfusionCase(dataDir, resDir, patientID, study_date)
% PerfCases = PerformGadgetronRecon_FindLGE_PerfusionCase('I:\BARTS', 'I:\ReconResults\BARTS', patientID, study_date)
% given a patientID and study date, find perf cases

studyDir = fullfile(resDir, study_date);

[subdirs, num] = FindSubDirs(studyDir);

PerfCases = [];

for ii=1:num
    
    v1 = strfind(subdirs{ii}, 'Perfusion');
    v2 = strfind(subdirs{ii}, patientID);
    v3 = strfind(subdirs{ii}, 'dicom');
    v4 = strfind(subdirs{ii}, 'Figure');
    
    if(~isempty(v1) & ~isempty(v2) & isempty(v3) & isempty(v4))
         PerfCases = [PerfCases; {subdirs{ii}}];      
    end    
end
