function copy_training_data_ai_cine_training_retro_binning_rt(aiDir, dstDir, pt_ids, files_record_picked, record)
% copy_training_data_ai_cine_training_retro_binning_rt(aiDir, dstDir, pt_ids, files_record_picked, record)

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    
    if(~isempty(record) && numel(pt_ids)==size(record.case_4chs, 1))
        case_4ch = record.case_4chs{pt};
        case_2ch = record.case_2chs{pt};
        case_3ch = record.case_3chs{pt};
        case_sax = record.case_saxs{pt};
        case_lvot = record.case_lvots{pt};
        case_aov = record.case_aovs{pt};
        case_rv = record.case_rvs{pt};
        case_aorticarch = record.case_aorticarchs{pt};
    else
        case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '4ch');
        case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '2ch');
        case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '3ch');
        case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'sa');
        case_lvot = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'lvot');
        case_aov = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'aov');
        case_rv = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'rv');
        case_aorticarch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'aorticarch');
    end

    try
        study_date = case_4ch.study_dates(end,:);
    catch
        try
            study_date = case_2ch.study_dates(end,:);
        catch
            try
                study_date = case_sax.study_dates(end,:);
            catch
                try
                    study_date = case_3ch.study_dates(end,:);
                catch
                    try
                        study_date = case_lvot.study_dates(end,:);
                    catch
                        try
                            study_date = case_aov.study_dates(end,:);
                        catch
                            study_date = case_rv.study_dates(end,:);
                        end
                    end
                end
            end
        end
    end
         
    if(size(case_4ch, 1)>0) 
        process_one_view_all(aiDir, dstDir, case_4ch, 'ch4');
    end

    % -----------------------------

    if(size(case_3ch, 1)>0)
        process_one_view_all(aiDir, dstDir, case_3ch, 'ch3');
    end
    % -----------------------------

    if(size(case_2ch, 1)>0)
        process_one_view_all(aiDir, dstDir, case_2ch, 'ch2');
    end       
end
end

function process_one_view_all(dataDir, dstDir, case_4ch, view_str)
    for ii=1:size(case_4ch, 1)
        try
            fname = case_4ch.file_names{ii};
            case_4ch_dir = fullfile(dataDir, case_4ch.study_dates(ii,:), case_4ch.file_names{ii});
            [path, sname, ext] = fileparts(case_4ch_dir); 
            sname= sname(~isspace(sname));
            case_prefix = [view_str '_' sname];
            [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(fname);

            src_dir = fullfile(dataDir, case_4ch.study_dates(ii,:), patientID, case_prefix);

            dst_dir_4ch = fullfile(dstDir, case_4ch.study_dates(ii,:), patientID, case_prefix);
            mkdir(dst_dir_4ch);

            copyfile(fullfile(src_dir, '*.jpg'), dst_dir_4ch);
            copyfile(fullfile(src_dir, 'data.npy'), dst_dir_4ch);
            copyfile(fullfile(src_dir, 'data_acq.npy'), dst_dir_4ch);
            copyfile(fullfile(src_dir, 'acq_time.npy'), dst_dir_4ch);
            copyfile(fullfile(src_dir, 'physio_time.npy'), dst_dir_4ch);
            copyfile(fullfile(src_dir, 'record_header.mat'), dst_dir_4ch);
        catch
            rmdir(dst_dir_4ch);
        end
    end
end
