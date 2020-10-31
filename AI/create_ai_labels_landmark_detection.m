function create_ai_labels_landmark_detection(aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, prefix, suffix, use_3D, local_pic_dir)
% create_ai_labels_landmark_detection(aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, prefix, suffix, use_3D, local_pic_dir)

    if(~isempty(local_pic_dir))
        pic_dir = fullfile(aiDir, local_pic_dir);
    else
        pic_dir = fullfile(aiDir, 'jpg_pics');
    end
    
    ch4_pic_dir = fullfile(pic_dir, 'ch4')
    ch2_pic_dir = fullfile(pic_dir, 'ch2')
    ch3_pic_dir = fullfile(pic_dir, 'ch3')
    sax_pic_dir = fullfile(pic_dir, 'sax')
    
    ch2_label_file = fullfile(pic_dir, [prefix '_ch2_label_file.csv']);
    ch3_label_file = fullfile(pic_dir, [prefix '_ch3_label_file.csv']);
    ch4_label_file = fullfile(pic_dir, [prefix '_ch4_label_file.csv']);
    sax_label_file = fullfile(pic_dir, [prefix '_sax_label_file.csv']);
    
    phases = [1:10:30, 13, 30];

    pts_loaded_ch4 = [];
    pts_loaded_ch2 = [];
    pts_loaded_ch3 = [];
    ind_ch4 = 1;
    ind_ch2 = 1;
    ind_ch3 = 1;
    
    for pt=1:numel(pt_ids)

        closeall

        pt_id = pt_ids{pt};

        disp(['-------> ' num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);

        if(numel(pt_ids)==size(case_4chs, 1))
            case_4ch = case_4chs{pt};
            case_2ch = case_2chs{pt};
            case_3ch = case_3chs{pt};
            case_sax = case_saxs{pt};
        else
            case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, {'4ch', 'hla'});
            case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, {'2ch', 'vla'});
            case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '3ch');
        end
        
        if(numel(case_4ch)==0 & numel(case_2ch)==0 & numel(case_3ch)==0)
            continue;
        end

        try
            if(numel(case_4ch)>0)
                study_date = case_4ch.study_dates(end,:);
            end
            if(numel(case_2ch)>0)
                study_date = case_2ch.study_dates(end,:);
            end
            if(numel(case_3ch)>0)
                study_date = case_3ch.study_dates(end,:);
            end

            case_prefix = [study_date '_' pt_id]
        catch
            case_prefix = [pt_id]
        end
        
        dst_dir = fullfile(aiDir, pt_id);
        if(~exist(dst_dir))
            dst_dir = fullfile(aiDir, study_date, pt_id);
        end
        
        for d=1:size(case_4ch,1)            
            case_4ch_dir = fullfile(aiDir, case_4ch.file_names{d});
            [path, sname, ext] = fileparts(case_4ch_dir);
            sname= sname(~isspace(sname));
            dst_dir_4ch = fullfile(dst_dir, ['ch4_' sname]);
            
            if(use_3D)
                ch4_pts_file = fullfile(dst_dir, 'res_ai', ['CH4_AI_pts_3D_' sname '_' suffix '.npy']);
            else
                ch4_pts_file = fullfile(dst_dir, 'res_ai', ['CH4_AI_pts_' sname '_' suffix '.npy']);
            end
            pts_loaded_ch4{ind_ch4, 1} = str2num(pt_id);
            pts_loaded_ch4{ind_ch4, 2} = sname;
            
            if(exist(ch4_pts_file))
                disp(['load ' ch4_pts_file]);
                pts_loaded_ch4{ind_ch4, 3} = readNPY(ch4_pts_file);
            else
                pts_loaded_ch4{ind_ch4, 3} = [];
            end
            
            ind_ch4 = ind_ch4 + 1;
        end  
        
        for d=1:size(case_3ch,1)            
            case_3ch_dir = fullfile(aiDir, case_3ch.file_names{d});
            [path, sname, ext] = fileparts(case_3ch_dir);
            sname= sname(~isspace(sname));
            dst_dir_3ch = fullfile(dst_dir, ['ch3_' sname]);
            
            if(use_3D)
                ch3_pts_file = fullfile(dst_dir, 'res_ai', ['CH3_AI_pts_3D_' sname '_' suffix '.npy']);
            else
                ch3_pts_file = fullfile(dst_dir, 'res_ai', ['CH3_AI_pts_' sname '_' suffix '.npy']);
            end
            pts_loaded_ch3{ind_ch3, 1} = str2num(pt_id);
            pts_loaded_ch3{ind_ch3, 2} = sname;
            
            if(exist(ch3_pts_file))
                disp(['load ' ch3_pts_file]);
                pts_loaded_ch3{ind_ch3, 3} = readNPY(ch3_pts_file);
            else
                pts_loaded_ch3{ind_ch3, 3} = [];
            end
            
            ind_ch3 = ind_ch3 + 1;
        end 
        
        for d=1:size(case_2ch,1)            
            case_2ch_dir = fullfile(aiDir, case_2ch.file_names{d});
            [path, sname, ext] = fileparts(case_2ch_dir);
            sname= sname(~isspace(sname));
            dst_dir_2ch = fullfile(dst_dir, ['ch2_' sname]);
            
            if(use_3D)
                ch2_pts_file = fullfile(dst_dir, 'res_ai', ['CH2_AI_pts_3D_' sname '_' suffix '.npy']);
            else
                ch2_pts_file = fullfile(dst_dir, 'res_ai', ['CH2_AI_pts_' sname '_' suffix '.npy']);
            end
            
            pts_loaded_ch2{ind_ch2, 1} = str2num(pt_id);
            pts_loaded_ch2{ind_ch2, 2} = sname;
            
            if(exist(ch2_pts_file))
                disp(['load ' ch2_pts_file]);
                pts_loaded_ch2{ind_ch2, 3} = readNPY(ch2_pts_file);
            else
                pts_loaded_ch2{ind_ch2, 3} = [];
            end
            
            ind_ch2 = ind_ch2 + 1;
        end 
    end
    
    process_one_view(ch2_pic_dir, ch2_label_file, pts_loaded_ch2(:,1), pts_loaded_ch2(:,2), pts_loaded_ch2(:,3));
    if(size(pts_loaded_ch3,1)>0)
        process_one_view(ch3_pic_dir, ch3_label_file, pts_loaded_ch3(:,1), pts_loaded_ch3(:,2), pts_loaded_ch3(:,3));
    end
    if(size(pts_loaded_ch4,1)>0)
        process_one_view(ch4_pic_dir, ch4_label_file, pts_loaded_ch4(:,1), pts_loaded_ch4(:,2), pts_loaded_ch4(:,3));
    end
    
    fclose('all')
end

function process_one_view(pic_dir, label_file, pt_ids, snames, pts)
%     filename,file_size,file_attributes,region_count,region_id,region_shape_attributes,region_attributes
%     20100722_1408_ch2_ED.jpg,11142,"{""Breathing Artifacts"":""no"",""Blurred myocardium"":""no"",""Low SNR"":""no"",""Wrong CMR View"":""no""}",3,0,"{""name"":""point"",""cx"":162,""cy"":107}","{}"
%     20100722_1408_ch2_ED.jpg,11142,"{""Breathing Artifacts"":""no"",""Blurred myocardium"":""no"",""Low SNR"":""no"",""Wrong CMR View"":""no""}",3,1,"{""name"":""point"",""cx"":142,""cy"":133}","{}"

    fid = fopen(label_file, 'w+');
    try
        fprintf(fid, ['filename,file_size,file_attributes,region_count,region_id,region_shape_attributes,region_attributes\n']);
%         if(ishandle(fid))
            [fnames, num] = findFILE(pic_dir, '*.jpg');
            for n=1:num
                [path, fname, ext] = fileparts(fnames{n});

                if(mod(n, 1000)==0)
                    fname
                end
                
                ind = find(fname=='_');

                pt_ind_ind = -1;
                
                if(numel(snames)>0)

                    try
                        ind2 = find(fname=='-');
                        fname_ismrmrd = fname(ind(1)+1:ind2(1)+6);
                        [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(fname_ismrmrd);
                        pt_id = str2num(patientID);
                         
                        if(~isempty(strfind(fname, 'ED')))
                            phs = 1;
                        else
                            phs = 13;
                        end
                    catch
                        pt_id = fname(1:ind(1)-1);
                        pt_id = str2num(pt_id);
                        phs = fname(ind(end)+1:end);
                        phs = str2num(phs);
                    end
                    
                    for pp=1:size(pt_ids,1)
                        if(pt_ids{pp}==pt_id & ~isempty(strfind(fnames{n}, snames{pp})))
                            pt_ind_ind = pp;
                            break;
                        end
                    end
                else
                    
                    study_date = fname(1:ind(1)-1);
                    pt_id = fname(ind(1)+1:ind(2)-1);
                    pt_id = str2num(pt_id);
                    phs = fname(ind(end)+1:end);
                    phs = str2num(phs);
                
                    for pp=1:size(pt_ids,1)
                        if(pt_ids{pp}==pt_id)
                            pt_ind_ind = pp;
                            break;
                        end
                    end
                end
                
                info = dir(fnames{n});
                if(pt_ind_ind<0)
                    content_str = [fname ext ',' num2str(info.bytes) ',"{""Breathing Artifacts"":""no"",""Blurred myocardium"":""no"",""Low SNR"":""no"",""Wrong CMR View"":""no""}"'];                
                    pt_str = [content_str ',' num2str(0) ',' num2str(0) ', "{}","{}"'];
                    fprintf(fid, [pt_str '\n']);               
                    continue;
                end
                
                curr_pt = pts{pt_ind_ind};
                if(isempty(curr_pt))                   
                    content_str = [fname ext ',' num2str(info.bytes) ',"{""Breathing Artifacts"":""no"",""Blurred myocardium"":""no"",""Low SNR"":""no"",""Wrong CMR View"":""no""}"'];                
                    pt_str = [content_str ',' num2str(0) ',' num2str(0) ', "{}","{}"'];
                    fprintf(fid, [pt_str '\n']);               
                    continue;
                end
                curr_pt = curr_pt(:,:,phs);
                
                N_pt = size(curr_pt, 1);
                num_pt = 0;
                for p=1:N_pt                    
                    if(curr_pt(p,1)>0)
                        num_pt = num_pt + 1;
                    end
                end
                
                content_str = [fname ext ',' num2str(info.bytes) ',"{""Breathing Artifacts"":""no"",""Blurred myocardium"":""no"",""Low SNR"":""no"",""Wrong CMR View"":""no""}"'];
                
                for p=1:N_pt      
                    if(curr_pt(p,1)>0)
                        pt_str = [content_str ',' num2str(num_pt) ',' num2str(p-1) ', "{""name"":""point"",""cx"":' num2str(round(curr_pt(p,1))) ',""cy"":' num2str(round(curr_pt(p,2))) '}","{}"'];
                        fprintf(fid, [pt_str '\n']);
                    end
                end
            end
%         end
        fclose(fid);
    catch e
        fclose(fid);
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
    end
end
