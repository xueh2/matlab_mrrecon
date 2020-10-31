
function files_record_selected = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record, patient_id, prot, only_one_case_needed)
% files_record_selected = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record, patient_id, prot)
% files_record_selected = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record, '32616742', '4ch')

if(nargin<4)
    only_one_case_needed = 0;
end

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
    %if(strcmp(patient_id, files_record.patientIDs{ii})==1)
    if(isnumeric(patient_id))
        has_pt = (patient_id==str2num(files_record.patientIDs{ii}));
    else
        has_pt = (strcmp(patient_id, files_record.patientIDs{ii}));
    end
    if(has_pt)
        if(iscell(prot))
            for k=1:numel(prot)
                try
                    if(~isempty(strfind(lower(files_record.prots{ii}), lower(prot{k}))) | ~isempty(strfind(lower(files_record.view_strs(ii, :)), lower(prot{k}))))
                        ind = [ind; ii];
                        if(only_one_case_needed)
                            break;
                        end
                    end
                catch
                    if(~isempty(strfind(lower(files_record.prots{ii}), lower(prot{k}))))
                        ind = [ind; ii];
                        if(only_one_case_needed)
                            break;
                        end
                    end
                end
            end
        else
            try
                if(~isempty(strfind(lower(files_record.prots{ii}), lower(prot))) | ~isempty(strfind(lower(files_record.view_strs(ii, :)), lower(prot))))
                    ind = [ind; ii];
                    if(only_one_case_needed)
                        break;
                    end
                end
            catch
                if(~isempty(strfind(lower(files_record.prots{ii}), lower(prot))))
                    ind = [ind; ii];
                    if(only_one_case_needed)
                        break;
                    end
                end
            end
        end
    end        
end

if(numel(ind)>0)
    if(only_one_case_needed)
        files_record_selected = files_record(ind(1), :);
    else
        files_record_selected = files_record(ind, :);
    end
end
