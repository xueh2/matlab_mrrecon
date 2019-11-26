
function files_record_selected = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record, patient_id, prot)
% files_record_selected = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record, patient_id, prot)
% files_record_selected = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record, '32616742', '4ch')


file_names = [];
study_dates = [];
study_times = [];
patientIDs = [];
studyIDs = [];
scannerIDs = [];
prots = [];
headers = [];
measurementIDs = [];

ind = [];
files_record_selected = [];

% check protocols
for ii=1:size(files_record,1)                
    if(strcmp(patient_id, files_record.patientIDs{ii})==1)
        if(~isempty(strfind(lower(files_record.prots{ii}), lower(prot))))
            ind = [ind; ii];
        end
    end        
end

if(numel(ind)>0)
    files_record_selected = files_record(ind, :);
end
