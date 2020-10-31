function create_ai_labels_landmark_detection_Perf(aiDir, pt_ids, files_record_picked, prefix, suffix, check_processed)
% create_ai_labels_landmark_detection_Perf(aiDir, pt_ids, files_record_picked, prefix, suffix, check_processed)

    pic_dir = fullfile(aiDir, 'jpg_pics');
    stress_pic_dir = fullfile(pic_dir, 'stress_eigen_image')
    rest_pic_dir = fullfile(pic_dir, 'rest_eigen_image')
    stress_pd_pic_dir = fullfile(pic_dir, 'stress_pd_image')
    rest_pd_pic_dir = fullfile(pic_dir, 'rest_pd_image')
    
    stress_label_file = fullfile(pic_dir, [prefix '_stress_label_file.csv']);
    rest_label_file = fullfile(pic_dir, [prefix '_rest_label_file.csv']);
    stress_pd_label_file = fullfile(pic_dir, [prefix '_stress_pd_label_file.csv']);
    rest_pd_label_file = fullfile(pic_dir, [prefix '_rest_pd_label_file.csv']);
    
    pts_loaded = cell(numel(pt_ids), 4);
    
    for pt=1:numel(pt_ids)

        closeall

        pt_id = pt_ids{pt};

        disp(['-------> ' num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);

        case_stress = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'stress', 1);
        case_rest = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'rest', 1);
        
        if(numel(case_stress)==0 & numel(case_rest)==0)
            continue;
        end

        if(numel(case_stress)>0)
            study_date = case_stress.study_dates(end,:);
        end
        if(numel(case_rest)>0)
            study_date = case_rest.study_dates(end,:);
        end
        
        case_prefix = [study_date '_' pt_id]
        dst_dir = fullfile(aiDir, study_date, pt_id);
        
        stress_pts_file = fullfile(dst_dir, 'res_ai', ['stress_AI_pts' suffix '.npy']);
        stress_pd_pts_file = fullfile(dst_dir, 'res_ai', ['stress_pd_AI_pts' suffix '.npy']);
        
        rest_pts_file = fullfile(dst_dir, 'res_ai', ['rest_AI_pts' suffix '.npy']);
        rest_pd_pts_file = fullfile(dst_dir, 'res_ai', ['rest_pd_AI_pts' suffix '.npy']);
        
%         if(~exist(stress_pts_file) | ~exist(rest_pts_file) | ~exist(stress_pd_pts_file) | ~exist(rest_pd_pts_file))
%             continue;
%         end
        
        pts_loaded{pt, 1} = str2num(pt_id);
        
        if(exist(stress_pts_file))
            pts_loaded{pt, 2} = readNPY(stress_pts_file);
        end
        
        if(exist(rest_pts_file))
            pts_loaded{pt, 3} = readNPY(rest_pts_file);
        end
        
        if(exist(stress_pd_pts_file))
            a = readNPY(stress_pd_pts_file);
            pts_loaded{pt, 4} = a;
        end
        
        if(exist(rest_pd_pts_file))
            a = readNPY(rest_pd_pts_file);
            pts_loaded{pt, 5} = a;
        end
    end
    
    process_one_view(stress_pic_dir, stress_label_file, pts_loaded(:,1), pts_loaded(:,2));
    process_one_view(rest_pic_dir, rest_label_file, pts_loaded(:,1), pts_loaded(:,3));
    
    process_one_view(stress_pd_pic_dir, stress_pd_label_file, pts_loaded(:,1), pts_loaded(:,4));
    process_one_view(rest_pd_pic_dir, rest_pd_label_file, pts_loaded(:,1), pts_loaded(:,5));
   
    fclose('all')
end

function process_one_view(pic_dir, label_file, pt_ids, pts)
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
                study_date = fname(1:ind(1)-1);
                pt_id = fname(ind(1)+1:ind(2)-1);
                pt_id = str2num(pt_id);
                phs = fname(ind(end)+1:end);
                phs = str2num(phs);

                pt_ind_ind = -1;
                for pp=1:size(pt_ids,1)
                    if(pt_ids{pp}==pt_id)
                        pt_ind_ind = pp;
                        break;
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
    catch
        fclose(fid);
    end
end
