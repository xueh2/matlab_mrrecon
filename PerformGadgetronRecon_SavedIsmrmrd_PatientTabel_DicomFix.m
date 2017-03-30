
function PerformGadgetronRecon_SavedIsmrmrd_PatientTabel_DicomFix(resDir, patient_table, dat_lines, dat_columns)
% PerformGadgetronRecon_SavedIsmrmrd_PatientTabel_DicomFix(dataDir, patient_table, dat_lines, dat_columns)
% setenv('OutputFormat', 'h5')

getenv('GT_HOST')
getenv('GT_PORT')

GTHome = getenv('GADGETRON_HOME')
GTConfigFolder = fullfile(GTHome, 'share/gadgetron/config');
date_suffix = datestr(date, 'yyyymmdd');

xmlUsed = '%GADGETRON_DIR%\install\schema/IsmrmrdParameterMap_Siemens_Perfusion.xml';

% ------------------------------------------------------------

% find data

files = [];

configNames = [];
study_dates = [];
study_times = [];

num = numel(dat_lines);
num_dat = numel(dat_columns);

for n=1:num  
    for kk=1:num_dat
        
        [pathstr, name, ext] = fileparts(patient_table{dat_lines(n), dat_columns(kk)});

        % find scanner ID, patient ID, study ID, measurement ID, study date and time
        [configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name);
        
        if( str2num(measurementID) > 10000 )
            continue;
        end

        dicomDir = fullfile(resDir, study_date, [name '_dicom']);
        
        str_to_find = {'_Vascular_Volume_Map', 'PS_Map', 'Perf_Map', 'Flow_Map_SLC', 'Gd_Extraction_Map', 'Interstitial_Volume_Map'};
        map = PerfColorMap;

        for kk=1:numel(str_to_find)
        %     [Vb, numVb] = findFILE(dstDir, '*_Vascular_Volume_Map*');
        %     [PS, numPS] = findFILE(dstDir, '*PS_Map*');
        %     [Ki, numKi] = findFILE(dstDir, '*Perf_Map*');
        %     [F, numF] = findFILE(dstDir, '*Flow_Map*');
        %     [E, numE] = findFILE(dstDir, '*Gd_Extraction_Map*');
        %     [Visf, numVisf] = findFILE(dstDir, '*Interstitial_Volume_Map*');

            [Vb, numVb] = findFILE(dicomDir, ['*' str_to_find{kk} '*'])

            for n=1:numVb    
                info = dicominfo(Vb{n});
                x = dicomread(Vb{n});    
                dicomwrite(x, map, Vb{n}, info);    
            end
        end
    end
end 
