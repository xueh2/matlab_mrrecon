function [record_4CH, record_3CH, record_2CH] = compare_Barts_Cine_GLS_results(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, suffix)
% record = compare_Barts_Cine_GLS_results(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, suffix)

    pt_ids_4CH = [];
    snames_4CH = [];
    pt_ids_3CH = [];
    snames_3CH = [];
    pt_ids_2CH = [];
    snames_2CH = [];
    
    GLS_4CH = [];
    GLS_3CH = [];
    GLS_2CH = [];
    GLS_4CH_3D = [];
    GLS_3CH_3D = [];
    GLS_2CH_3D = [];
    
    MAPSE_4CH = [];
    MAPSE_3CH = [];
    MAPSE_2CH = [];
    MAPSE_4CH_3D = [];
    MAPSE_3CH_3D = [];
    MAPSE_2CH_3D = [];
    
    LV_length_PT1_4CH = [];
    LV_length_PT1_3CH = [];
    LV_length_PT1_2CH = [];
    LV_length_PT1_4CH_3D = [];
    LV_length_PT1_3CH_3D = [];
    LV_length_PT1_2CH_3D = [];

    LV_length_PT2_4CH = [];
    LV_length_PT2_3CH = [];
    LV_length_PT2_2CH = [];
    LV_length_PT2_4CH_3D = [];
    LV_length_PT2_3CH_3D = [];
    LV_length_PT2_2CH_3D = [];
    
    pts_4CH = [];
    pts_3CH = [];
    pts_2CH = [];
    pts_4CH_3D = [];
    pts_3CH_3D = [];
    pts_2CH_3D = [];
    
    for pt=1:numel(pt_ids)

        closeall

        try
            pt_id = pt_ids{pt};
        catch
            pt_id = pt_ids(pt);
        end
        if(isnumeric(pt_id))
            pt_id = num2str(pt_id);
        end
        
        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);

        if(numel(pt_ids)==numel(case_4chs))
            case_4ch = case_4chs{pt};
            case_2ch = case_2chs{pt};
            case_3ch = case_3chs{pt};
            case_sax = case_saxs{pt};
        else
            case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '4ch', 0);
            case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '2ch', 0);
            case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '3ch', 0);
            case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'sa', 0);
            
            if(numel(case_4ch)==0 & numel(case_2ch)==0 & numel(case_3ch)==0)
                case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'CH4', 0);
                case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'CH2', 0);
                case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'CH3', 0);
                case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'SAX', 0);
            end
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
            dst_dir = fullfile(aiDir, study_date, pt_id);
        catch
            case_prefix = [pt_id]
            dst_dir = fullfile(aiDir, pt_id);
        end
        
        if(~exist(dst_dir))
            dst_dir = fullfile(aiDir, pt_id);
        end
        
        % 4ch
        dst_dir_4ch = fullfile(dst_dir, 'ch4');

        % 2ch
        dst_dir_2ch = fullfile(dst_dir, 'ch2');

        % 3ch
        dst_dir_3ch = fullfile(dst_dir, 'ch3');

        %sax
        dst_dir_sax = fullfile(dst_dir, 'sax');

        contourDir = fullfile(dst_dir, 'res_ai');

        ES_as_first_phase = 0;
        
        for d=1:size(case_4ch,1)
            case_4ch_dir = fullfile(resDir, case_4ch.file_names{d});
            [path, sname, ext] = fileparts(case_4ch_dir);
            sname= sname(~isspace(sname));
            dst_dir_4ch = fullfile(dst_dir, ['ch4_' sname]);
            
            pts_file = fullfile(contourDir, ['CH4_AI_pts' '_' sname '_' suffix '.npy']);
            
            pt_ids_4CH = [pt_ids_4CH; str2num(pt_id)];
            snames_4CH = [snames_4CH; {sname}];

            GLS = -1;
            GLS_3D = -1;
            MAPSE = -1;
            MAPSE_3D = -1;
            
            if(exist(pts_file))
                pts = readNPY(pts_file);                
                clear GLS MAPSE LV_Lengh_PT1 LV_Lengh_PT2
                [GLS, MAPSE, LV_Lengh_PT1, LV_Lengh_PT2] = compute_GLS(pts, ES_as_first_phase);
                save(fullfile(contourDir, ['CH4_AI_' sname '.mat']), 'GLS', 'MAPSE', 'LV_Lengh_PT1', 'LV_Lengh_PT2', 'pts');
            end
            
            LV_length_PT1_4CH = [LV_length_PT1_4CH; {LV_Lengh_PT1}];
            LV_length_PT2_4CH = [LV_length_PT2_4CH; {LV_Lengh_PT2}];   
            pts_4CH = [pts_4CH; {pts}];
            
            pts_file = fullfile(contourDir, ['CH4_AI_pts_3D' '_' sname '_' suffix '.npy']);
            if(exist(pts_file))
                pts = readNPY(pts_file);                
                clear GLS_3D MAPSE_3D LV_Lengh_PT1 LV_Lengh_PT2
                [GLS_3D, MAPSE_3D, LV_Lengh_PT1, LV_Lengh_PT2] = compute_GLS(pts, ES_as_first_phase);
                save(fullfile(contourDir, ['CH4_AI_3D_' sname '.mat']), 'GLS_3D', 'MAPSE_3D', 'LV_Lengh_PT1', 'LV_Lengh_PT2', 'pts');
            end
            
            GLS_4CH = [GLS_4CH;  max(GLS)];
            GLS_4CH_3D = [GLS_4CH_3D;  max(GLS_3D)];
            MAPSE_4CH = [MAPSE_4CH;  max(MAPSE)];
            MAPSE_4CH_3D = [MAPSE_4CH_3D;  max(MAPSE_3D)];
            
            LV_length_PT1_4CH_3D = [LV_length_PT1_4CH_3D; {LV_Lengh_PT1}];
            LV_length_PT2_4CH_3D = [LV_length_PT2_4CH_3D; {LV_Lengh_PT2}];            
            pts_4CH_3D = [pts_4CH_3D; {pts}];
        end
        
        for d=1:size(case_3ch,1)
            case_3ch_dir = fullfile(resDir, case_3ch.file_names{d});
            [path, sname, ext] = fileparts(case_3ch_dir);
            sname= sname(~isspace(sname));
            dst_dir_3ch = fullfile(dst_dir, ['ch3_' sname]);
            
            pts_file = fullfile(contourDir, ['CH3_AI_pts' '_' sname '_' suffix '.npy']);
            
            pt_ids_3CH = [pt_ids_3CH; str2num(pt_id)];
            snames_3CH = [snames_3CH; {sname}];

            GLS = -1;
            GLS_3D = -1;
            MAPSE = -1;
            MAPSE_3D = -1;
            
            if(exist(pts_file))
                pts = readNPY(pts_file);  
                clear GLS MAPSE LV_Lengh_PT1 LV_Lengh_PT2
                [GLS, MAPSE, LV_Lengh_PT1, LV_Lengh_PT2] = compute_GLS(pts, ES_as_first_phase);
                save(fullfile(contourDir, ['CH3_AI_' sname '.mat']), 'GLS', 'MAPSE', 'LV_Lengh_PT1', 'LV_Lengh_PT2', 'pts');
            end
            
            LV_length_PT1_3CH = [LV_length_PT1_3CH; {LV_Lengh_PT1}];
            LV_length_PT2_3CH = [LV_length_PT2_3CH; {LV_Lengh_PT2}];
            pts_3CH = [pts_3CH; {pts}];
            
            pts_file = fullfile(contourDir, ['CH3_AI_pts_3D' '_' sname '_' suffix '.npy']);
            if(exist(pts_file))
                pts = readNPY(pts_file);                
                clear GLS_3D MAPSE_3D LV_Lengh_PT1 LV_Lengh_PT2
                [GLS_3D, MAPSE_3D, LV_Lengh_PT1, LV_Lengh_PT2] = compute_GLS(pts, ES_as_first_phase);
                save(fullfile(contourDir, ['CH3_AI_3D_' sname '.mat']), 'GLS_3D', 'MAPSE_3D', 'LV_Lengh_PT1', 'LV_Lengh_PT2', 'pts');
            end
            
            GLS_3CH = [GLS_3CH;   max(GLS)];
            GLS_3CH_3D = [GLS_3CH_3D;   max(GLS_3D)];  
            MAPSE_3CH = [MAPSE_3CH;  max(MAPSE)];
            MAPSE_3CH_3D = [MAPSE_3CH_3D;  max(MAPSE_3D)];
            
            LV_length_PT1_3CH_3D = [LV_length_PT1_3CH_3D; {LV_Lengh_PT1}];
            LV_length_PT2_3CH_3D = [LV_length_PT2_3CH_3D; {LV_Lengh_PT2}];
            pts_3CH_3D = [pts_3CH_3D; {pts}];
        end
        
        for d=1:size(case_2ch,1)
            case_2ch_dir = fullfile(resDir, case_2ch.file_names{d});
            [path, sname, ext] = fileparts(case_2ch_dir);
            sname= sname(~isspace(sname));
            dst_dir_2ch = fullfile(dst_dir, ['ch2_' sname]);
            
            pts_file = fullfile(contourDir, ['CH2_AI_pts' '_' sname '_' suffix '.npy']);
            
            pt_ids_2CH = [pt_ids_2CH; str2num(pt_id)];
            snames_2CH = [snames_2CH; {sname}];

            GLS = -1;
            GLS_3D = -1;
            MAPSE = -1;
            MAPSE_3D = -1;
            
            if(exist(pts_file))
                pts = readNPY(pts_file);                
                if(size(pts, 3)>10)
                    clear GLS MAPSE LV_Lengh_PT1 LV_Lengh_PT2
                    [GLS, MAPSE, LV_Lengh_PT1, LV_Lengh_PT2] = compute_GLS(pts, ES_as_first_phase);
                    save(fullfile(contourDir, ['CH2_AI_' sname '.mat']), 'GLS', 'MAPSE', 'LV_Lengh_PT1', 'LV_Lengh_PT2', 'pts');
                end
            end
                       
            LV_length_PT1_2CH = [LV_length_PT1_2CH; {LV_Lengh_PT1}];
            LV_length_PT2_2CH = [LV_length_PT2_2CH; {LV_Lengh_PT2}];
            
            pts_2CH = [pts_2CH; {pts}];
            
            pts_file = fullfile(contourDir, ['CH2_AI_pts_3D' '_' sname '_' suffix '.npy']);
            if(exist(pts_file))
                pts = readNPY(pts_file);                
                if(size(pts, 3)>10)
                    clear GLS_3D MAPSE_3D LV_Lengh_PT1 LV_Lengh_PT2
                    [GLS_3D, MAPSE_3D, LV_Lengh_PT1, LV_Lengh_PT2] = compute_GLS(pts, ES_as_first_phase);
                    save(fullfile(contourDir, ['CH2_AI_3D_' sname '.mat']), 'GLS_3D', 'MAPSE_3D', 'LV_Lengh_PT1', 'LV_Lengh_PT2', 'pts');
                end
            end
            
            GLS_2CH = [GLS_2CH;   max(GLS)];
            GLS_2CH_3D = [GLS_2CH_3D;   max(GLS_3D)];   
            MAPSE_2CH = [MAPSE_2CH;  max(MAPSE)];
            MAPSE_2CH_3D = [MAPSE_2CH_3D;  max(MAPSE_3D)];
            
            LV_length_PT1_2CH_3D = [LV_length_PT1_2CH_3D; {LV_Lengh_PT1}];
            LV_length_PT2_2CH_3D = [LV_length_PT2_2CH_3D; {LV_Lengh_PT2}];
            
            pts_2CH_3D = [pts_2CH_3D; {pts}];
        end        
    end
    
    record_4CH = table(pt_ids_4CH, snames_4CH, GLS_4CH, GLS_4CH_3D, MAPSE_4CH, MAPSE_4CH_3D, LV_length_PT1_4CH, LV_length_PT2_4CH, LV_length_PT1_4CH_3D, LV_length_PT2_4CH_3D, pts_4CH, pts_4CH_3D);
    record_3CH = table(pt_ids_3CH, snames_3CH, GLS_3CH, GLS_3CH_3D, MAPSE_3CH, MAPSE_3CH_3D, LV_length_PT1_3CH, LV_length_PT2_3CH, LV_length_PT1_3CH_3D, LV_length_PT2_3CH_3D, pts_3CH, pts_3CH_3D);
    record_2CH = table(pt_ids_2CH, snames_2CH, GLS_2CH, GLS_2CH_3D, MAPSE_2CH, MAPSE_2CH_3D, LV_length_PT1_2CH, LV_length_PT2_2CH, LV_length_PT1_2CH_3D, LV_length_PT2_2CH_3D, pts_2CH, pts_2CH_3D);
end

% function [GLS, MAPSE, LV_Lengh_PT1, LV_Lengh_PT2] = compute_GLS(pts, ES_as_first_phase)
% 
%     if(~isempty(find(pts<2)))
%         GLS = -1;
%         MAPSE = -1;
%         return
%     end
% 
%     GLS = -1;
%     MAPSE = -1;
%     
%     ptc = 0.5 * (pts(1,:,:) + pts(2,:,:));
%     ptc = squeeze(ptc);
%     
%     N = size(pts, 3);
%     
%     L = zeros(N,1);
%     LV_Lengh_PT1 = zeros(N, 1);
%     LV_Lengh_PT2 = zeros(N, 1);
%     
%     for j = 1:N        
%         L(j) = norm([ptc(1,j), ptc(2,j)] - [pts(3,1,j), pts(3,2,j)]);
%         LV_Lengh_PT1(j) = norm([pts(1,1,j), pts(1,2,j)] - [pts(3,1,j), pts(3,2,j)]);
%         LV_Lengh_PT2(j) = norm([pts(2,1,j), pts(2,2,j)] - [pts(3,1,j), pts(3,2,j)]);
%     end
% 
%     if(ES_as_first_phase)
%         GLS = 100 * (L(1) - L) / L(1);
%     else
%         maxL = max(L(:));
%         sorted_L = sort(L);
%         %GLS = 100 * (maxL - L) / mean(sorted_L(N-5:N));
%         % maxL = (0.5 * (L(1)+L(end)));
%         GLS = 100 * (maxL - L) / (0.5 * (L(1)+L(end)));
%     end
%     
%     MAPSE = zeros(N, 1);
%     ptc_ED = (ptc(:,1) + ptc(:,end)) / 2;
%     for j = 1:N        
%         MAPSE(j) = norm([ptc(1,j), ptc(2,j)] - [ptc_ED(1,1), ptc_ED(2,1)]);
%     end
% end
% 
