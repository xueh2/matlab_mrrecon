function collect_scan_info_ai(dataDir, resDir, aiDir, files_all)
% collect_scan_info_ai(dataDir, resDir, aiDir, files_all)

for n = 1:size(files_all,1)
    
    fname = files_all{n, 1};
    if(iscell(fname))
        fname = fname{1};
    end
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(fname);
    
    h5_name = fullfile(dataDir, study_dates, [fname '.h5']);
    
    case_dir = fullfile(resDir, study_dates, fname);   
    dst_dir = fullfile(aiDir, study_dates, fname);
    if ~exist(dst_dir)
        mkdir(dst_dir)
    end
                     
    try
        if (~exist(fullfile(dst_dir, 'ismrmrd_hdr.mat')))
            dset = ismrmrd.Dataset(h5_name, 'dataset');
            hdr = ismrmrd.xml.deserialize(dset.readxml);
            dset.close();

            disp(['--> ' num2str(n) ' out of ' num2str(size(files_all,1)) ' - ' num2str(hdr.acquisitionSystemInformation.systemFieldStrength_T) ' -- ' hdr.measurementInformation.protocolName ' -- ' hdr.acquisitionSystemInformation.institutionName '--' files_all{n}]);

            save(fullfile(dst_dir, 'ismrmrd_hdr.mat'), 'hdr');
        else
            disp(['--> already processed -- ' files_all{n}]);
        end
    catch
        disp(['--> failed for ' files_all{n}]);
    end
end