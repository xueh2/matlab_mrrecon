function [aif_scan_geometry_info, scan_geometry_info] = read_in_GT_Perf_DebugOutput_scan_geometry_info(resDir)

% read in Gadgetron perfusion debug output results
% [aif_scan_geometry_info, scan_geometry_info] = read_in_GT_Perf_DebugOutput_scan_geometry_info(resDir)

    try
        slc = 0;   
        for n=1:8

    %         filename = ['perf_moco_upsampled_' num2str(n-1) '.hdr'];
    %         if(~isFileExist(fullfile(resDir, 'DebugOutput', filename)))
    %             break;
    %         end        

            [names, num] = findFILE(resDir, ['*_SLC' num2str(n-1) '*103.attrib']);

            if(num==0)
                break;
            end
            slc = slc + 1;
        end
      
        MS_mode = 0;
        if(slc==1)
            slc = 0;
            for n=1:8
                [names, num] = findFILE(resDir, ['*_SLC' num2str(n-1) '*' num2str(n) '03.attrib']);

                if(num==0)
                    break;
                end
                slc = slc + 1;
                MS_mode = 1;
            end
        end

        disp(['Total ' num2str(slc) ' is found ...']);

        % read in aif slice and image position
        xmlContent = xml_load(fullfile(resDir, ['results_SLC0_CON0_PHS0_REP0_SET0_AVE0_1_1104.attrib']));
        aif_slice_dir = getXMLField(xmlContent, 'slice_dir');
        aif_read_dir = getXMLField(xmlContent, 'read_dir');
        aif_phase_dir = getXMLField(xmlContent, 'phase_dir');
        aif_PatientPosition = getXMLField(xmlContent, 'PatientPosition');
        aif_patient_table_position = getXMLField(xmlContent, 'patient_table_position');    

        for s=1:slc
            if(MS_mode==0)
                [names, num] = findFILE(resDir, ['*_SLC' num2str(s-1) '*103.attrib']);
            else
                [names, num] = findFILE(resDir, ['*_SLC' num2str(s-1) '*' num2str(s) '03.attrib']);
            end
            names{1}
            xmlContent = xml_load(names{1});
            slice_dir(s,:) = getXMLField(xmlContent, 'slice_dir');
            read_dir(s,:) = getXMLField(xmlContent, 'read_dir');
            phase_dir(s,:) = getXMLField(xmlContent, 'phase_dir');
            PatientPosition(s,:) = getXMLField(xmlContent, 'PatientPosition');
            patient_table_position(s,:) = getXMLField(xmlContent, 'patient_table_position');
        end

        aif_scan_geometry_info = table(aif_slice_dir, aif_read_dir, aif_phase_dir, aif_PatientPosition, aif_patient_table_position);    
        scan_geometry_info = table(slice_dir, read_dir, phase_dir, PatientPosition, patient_table_position);    
    catch
        disp(['Try h5 ... ']);
        
        [names, num] = findFILE(resDir, '*.h5');
        
        aif_headers = h5read(names{1}, '/results/image_1104/header');
        perf_headers = h5read(names{1}, '/results/image_124/header');
        
        slc = numel(perf_headers.image_index);
        
        aif_slice_dir = aif_headers.slice_dir(:,2)';
        aif_read_dir = aif_headers.read_dir(:,2)';
        aif_phase_dir = aif_headers.phase_dir(:,2)';
        aif_PatientPosition = aif_headers.position(:,2)';
        aif_patient_table_position = aif_headers.patient_table_position(:,2)';
        
        slice_dir = perf_headers.slice_dir';
        read_dir = perf_headers.read_dir';
        phase_dir = perf_headers.phase_dir';
        PatientPosition = perf_headers.position';
        patient_table_position = perf_headers.patient_table_position';

        aif_scan_geometry_info = table(aif_slice_dir, aif_read_dir, aif_phase_dir, aif_PatientPosition, aif_patient_table_position);
        scan_geometry_info = table(slice_dir, read_dir, phase_dir, PatientPosition, patient_table_position);    
    end

end

function v = getXMLField(xmlContent, vname)

    for ii=1:numel(xmlContent)        
        if(strcmp(xmlContent(ii).meta(1).name, vname)==1)            
            for j=1:numel(xmlContent(ii).meta)
                v(j) = str2double(xmlContent(ii).meta(j).value);
            end            
        end        
    end
end