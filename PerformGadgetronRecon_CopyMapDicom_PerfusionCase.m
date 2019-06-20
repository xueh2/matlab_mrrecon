
function missing_cases = PerformGadgetronRecon_CopyMapDicom_PerfusionCase(perf_cases, rest_cases, resDir, dicomDir)
% PerformGadgetronRecon_CopyMapDicom_PerfusionCase(perf_cases, rest_cases, resDir, dicomDir)

rsN = size(perf_cases, 1);
N = rsN*2 + size(rest_cases, 1);

cases = cell(N,1);

for k=1:rsN  
    cases{1+(k-1)*2} = perf_cases{k, 2};
    cases{2+(k-1)*2} = perf_cases{k, 3};    
end

for k=rsN*2+1:N   
    cases{k} = rest_cases{k-rsN*2};
end

mbv = MBVColorMap(0);
mbf = PerfColorMap(0);
ps = PSColorMap(0);

for n=1:3
    mbv_interp(:,n) = interp(mbv(:,n), 2^16/256);
    mbf_interp(:,n) = interp(mbf(:,n), 2^16/256);
    ps_interp(:,n) = interp(ps(:,n), 2^16/256);
end

% closeall
% for n=1:N
%            
%     stressCase = cases{n};
%     disp('=================================================');
%     [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(stressCase);
%     
%     stress_dicom = fullfile(resDir, study_dates, [stressCase '_dicom']);   
%     dst_stress = fullfile(dicomDir, study_dates, [stressCase '_dicom']);
%     dst = dst_stress;
%     
%     [fname, numVb] = findFILE(dst, ['*' 'Flow_Map_SLC' '*']);
%     
%     if(numVb==0)
%         continue;
%     end
%     
%     clear I
%     for kk=1:numVb
%         I(:,:,kk) = dicomread(fname{kk});
%     end
% 
%     disp([num2str(n) ' out of ' num2str(N) ' - ' dst]);
%     figure; imagescn( double(I)/100, [0 6], [1 numVb], [12]); PerfColorMap;
%     pause
%     closeall
% end
% 
% error('stop');

missing_cases = [];

for n=1:N
           
    stressCase = cases{n};
    disp('=================================================');
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(stressCase);
    
    stress_dicom = fullfile(resDir, study_dates, [stressCase '_dicom']);
    if(isempty(dicomDir))
        dst_stress = stress_dicom;
    else
        dst_stress = fullfile(dicomDir, study_dates, [stressCase '_dicom']);
        mkdir(dst_stress);
    end

    dst = dst_stress;
            
    % str_to_find = {'_Vascular_Volume_Map', 'PS_Map', 'Perf_Map', 'Flow_Map_SLC', 'Gd_Extraction_Map', 'Interstitial_Volume_Map'};
    % str_to_find = {'Perf_Map', 'Flow_Map_SLC', 'FIG'};
    % str_to_find = {'_Vascular_Volume_Map', 'Flow_Map_SLC'};
    % str_to_find = {'Flow_Map_SLC'};
    str_to_find = {'_Vascular_Volume_Map', 'Flow_Map_SLC', 'PS_Map', };
        
    for kk=1:numel(str_to_find)

        if(kk==1 & numel(str_to_find)>1)
            c = mbv_interp;
        elseif(kk==2)
            c = mbf_interp;
        else
            c = ps_interp;
        end
        
        [fname, numVb] = findFILE(stress_dicom, ['*' str_to_find{kk} '*']);
        
        if(numVb==0)
            missing_cases = [missing_cases; {stress_dicom}];
            continue;
        end
        
        has_fixed = 0;
        for ii=1:numVb
            [pathstr, name, ext] = fileparts(fname{ii});
            if(~isempty(strfind(name, 'fixed')))
                has_fixed = 1;
                break;
            end
        end

        uid = dicomuid;
        
        for ii=1:numVb    
            [pathstr, name, ext] = fileparts(fname{ii});
            if(has_fixed)
                if(isempty(strfind(name, '_colored')))
                    continue;
                end
            end

            if(~isempty(strfind(name, '_colored')))
                continue;
            end
                
            dst_fname = fullfile(dst, [name '.dcm']);

            s = dicomread(fname{ii});
            info = dicominfo(fname{ii});

            info.ColorType = 'grayscale';
            info.PhotometricInterpretation = 'MONOCHROME2';
            % dicomwrite(uint16(s), dst_fname, info);

            info.SeriesInstanceUID = uid;
            info.SeriesNumber = info.SeriesNumber*10;
%             info.ColorType = 'indexed';
%             info.PhotometricInterpretation = 'PALETTE COLOR';

            cc = 65535 * c;
            cc = uint16(cc);

            info.RedPaletteColorLookupTableData = cc(:,1);
            info.GreenPaletteColorLookupTableData = cc(:,2);
            info.BluePaletteColorLookupTableData = cc(:,3);
            
            info.RedPaletteColorLookupTableDescriptor = [65535; 0; 16];
            info.BluePaletteColorLookupTableDescriptor = [65535; 0; 16];
            info.GreenPaletteColorLookupTableDescriptor = [65535; 0; 16];
            info.RescaleIntercept=0;
            info.RescaleSlope=1;
            info.RescaleType='US';
            
            % dicomwrite(s, c, fullfile(dst, [name '_colored.dcm']), info);
            % copyfile(fname{ii}, dst_fname);
            
            try
                dicomwrite(uint16(s), fullfile(dst, [name '_colored.dcm']), info, 'CreateMode','copy'); 
            catch
            end
        end
    end
end

missing_cases

