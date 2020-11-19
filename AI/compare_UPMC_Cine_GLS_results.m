function record = compare_UPMC_Cine_GLS_results(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, res, suffix)
% record = compare_UPMC_Cine_GLS_results(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, res, suffix)

    GLS_4CH = [];
    GLS_3CH = [];
    GLS_2CH = [];
    GLS_4CH_ED_first_phase = [];
    GLS_3CH_ED_first_phase = [];
    GLS_2CH_ED_first_phase = [];
    MAPSE_4CH = [];
    MAPSE_3CH = [];
    MAPSE_2CH = [];
    lMAPSE_4CH = [];
    
    GLS_4CH_2D = [];
    GLS_3CH_2D = [];
    GLS_2CH_2D = [];
    GLS_4CH_ED_first_phase_2D = [];
    GLS_3CH_ED_first_phase_2D = [];
    GLS_2CH_ED_first_phase_2D = [];
    MAPSE_4CH_2D = [];
    MAPSE_3CH_2D = [];
    MAPSE_2CH_2D = [];
    lMAPSE_4CH_2D = [];
    
    GLS_Strain = [];
    
    EF = [];
    ESV = [];
    EDV = [];
    MASS = [];
    Infarction = [];
    NonIschemicScar = [];
    LGE = [];
    BSA = [];
    BMI = [];
    Height = [];
    Weight = [];
    Gender = [];
    Age = [];
    FELKER = [];
    patients = [];
    pt_dirs = [];
    pids = [];
    
    for pt=1:numel(pt_ids)

        pt_id = pt_ids{pt};

        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);

        % find pt_id
        has_pt_id = 0;
        for k=1:size(res,1)
            res_pt_id = res(k, 15);
            if(res_pt_id == str2num(pt_id))
                has_pt_id = 1;
                break;
            end
        end
        
        if(~has_pt_id)
            disp(['-----> cannot find pd_id = ' pt_id]);
            continue;
        end
        
        if(isnan(res(k, 13)))
            disp(['-----> cannot find valid strain GLS, ' pt_id]);
            continue;
        end
        
        res_record = res(k, :);
            
        if(isnumeric(pt_id))
            disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' num2str(pt_id)]);
        else
            disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
        end

        if(numel(pt_ids)==size(case_4chs, 1))
            case_4ch = case_4chs{pt};
            case_2ch = case_2chs{pt};
            case_3ch = case_3chs{pt};
            case_sax = case_saxs{pt};
        else
            case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '4ch', 1);
            case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '2ch', 1);
            case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, '3ch', 1);
            case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'sa', 1);
            
            if(numel(case_4ch)==0 & numel(case_2ch)==0 & numel(case_3ch)==0)
                case_4ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'CH4', 1);
                case_2ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'CH2', 1);
                case_3ch = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'CH3', 1);
                case_sax = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'SAX', 1);
            end
        end

        if(numel(case_4ch)==0 & numel(case_2ch)==0 & numel(case_3ch)==0)
            continue;
        end

        if(numel(case_4ch)>0)
            study_date = case_4ch.study_dates(end,:);
            pt_id = case_4ch.patientIDs{end,:};
        end
        if(numel(case_2ch)>0)
            study_date = case_2ch.study_dates(end,:);
            pt_id = case_2ch.patientIDs{end,:};
        end
        if(numel(case_3ch)>0)
            study_date = case_3ch.study_dates(end,:);
            pt_id = case_3ch.patientIDs{end,:};
        end

        dst_dir = fullfile(aiDir, study_date, pt_id);
        if(~exist(dst_dir))
            dst_dir = fullfile(aiDir, pt_id);
        end
        
        [case_dirs, num] = FindSubDirs(dst_dir);
        
        % 4ch
        dst_dir_4ch = [];
        for tt=1:num
            if(~isempty(strfind(case_dirs{tt}, 'ch4')))
                dst_dir_4ch = fullfile(dst_dir, case_dirs{tt});
            end
        end
        
        suffix_ch4 = suffix;
        if(~isempty(dst_dir_4ch))
            [path, sname, ext] = fileparts(dst_dir_4ch); 
            sname= sname(~isspace(sname));
            suffix_ch4 = [sname(4:end) '_' suffix];
        end
        
        % 2ch
        dst_dir_2ch = [];
        for tt=1:num
            if(~isempty(strfind(case_dirs{tt}, 'ch2')))
                dst_dir_2ch = fullfile(dst_dir, case_dirs{tt});
            end
        end

        suffix_ch2 = suffix;
        if(~isempty(dst_dir_2ch))
            [path, sname, ext] = fileparts(dst_dir_2ch); 
            sname= sname(~isspace(sname));
            suffix_ch2 = [sname(4:end) '_' suffix];
        end
        
        % 3ch
        dst_dir_3ch = [];
        for tt=1:num
            if(~isempty(strfind(case_dirs{tt}, 'ch3')))
                dst_dir_3ch = fullfile(dst_dir, case_dirs{tt});
            end
        end

        suffix_ch3 = suffix;
        if(~isempty(dst_dir_3ch))
            [path, sname, ext] = fileparts(dst_dir_3ch); 
            sname= sname(~isspace(sname));
            suffix_ch3 = [sname(4:end) '_' suffix];
        end
        
        %sax
        dst_dir_sax = fullfile(dst_dir, 'sax');

        contourDir = fullfile(dst_dir, 'res_ai');

        patients = [patients; str2num(pt_id)];
        pids = [pids; pt_id];
        
        GLS_Strain = [GLS_Strain; res(k, 13)];
        EF = [EF; res(k, 7)];
        ESV = [ESV; res(k, 9)];
        EDV = [EDV; res(k, 8)];
        MASS = [MASS; res(k, 10)];
        Infarction = [Infarction; res(k, 11)];
        NonIschemicScar = [NonIschemicScar; res(k, 12)];
        LGE = [LGE; res(k, 14)];
        BSA = [BSA; res(k, 6)];
        BMI = [BMI; res(k, 5)];
        Height = [Height; res(k, 4)];
        Weight = [Weight; res(k, 3)];
        Gender = [Gender; res(k, 2)];
        Age = [Age; res(k, 1)];
        FELKER = [FELKER; res(k, 16)];        
        pt_dirs = [pt_dirs; {contourDir}];
        
        ch4_pt_file = fullfile(contourDir, ['CH4_AI_pts_3D' suffix '.npy']);
        if(~exist(ch4_pt_file))
            ch4_pt_file = fullfile(contourDir, ['CH4_AI_pts_3D' suffix_ch4 '.npy']);
        end        
        
        if(exist(ch4_pt_file))
            pts = readNPY(ch4_pt_file);
            
            ES_as_first_phase = 1;
            [GLS, MAPSE, lMAPSE] = compute_GLS(pts, ES_as_first_phase);
            ES_as_first_phase = 0;
            [GLS2, MAPSE2, lMAPSE2] = compute_GLS(pts, ES_as_first_phase);
            
            GLS_4CH = [GLS_4CH; max(GLS2)];
            GLS_4CH_ED_first_phase = [GLS_4CH_ED_first_phase; max(GLS)];
            MAPSE_4CH = [MAPSE_4CH; max(MAPSE2)];
            lMAPSE_4CH = [lMAPSE_4CH; max(lMAPSE2)];
        else
            disp(['-----> cannot find CH4_AI_pts_3D, ' pt_id ' - ' study_date ' - ' suffix]);
            GLS_4CH = [GLS_4CH; -1];
            GLS_4CH_ED_first_phase = [GLS_4CH_ED_first_phase; -1];
            MAPSE_4CH = [MAPSE_4CH; -1];
        end
        
        ch2_pt_file = fullfile(contourDir, ['CH2_AI_pts_3D' suffix '.npy']);
        if(~exist(ch2_pt_file))
            ch2_pt_file = fullfile(contourDir, ['CH2_AI_pts_3D' suffix_ch2 '.npy']);
        end        
        
        if(exist(ch2_pt_file))
            pts = readNPY(ch2_pt_file);
            
            ES_as_first_phase = 1;
            [GLS, MAPSE, V] = compute_GLS(pts, ES_as_first_phase);
            ES_as_first_phase = 0;
            [GLS2, MAPSE2, V] = compute_GLS(pts, ES_as_first_phase);
            
            GLS_2CH = [GLS_2CH; max(GLS2)];
            GLS_2CH_ED_first_phase = [GLS_2CH_ED_first_phase; max(GLS)];            
            MAPSE_2CH = [MAPSE_2CH; max(MAPSE2)];
        else
            disp(['-----> cannot find CH2_AI_pts_3D, ' pt_id ' - ' study_date ' - ' suffix]);
            GLS_2CH = [GLS_2CH; -1];
            GLS_2CH_ED_first_phase = [GLS_2CH_ED_first_phase; -1];
            MAPSE_2CH = [MAPSE_2CH; -1];
        end
           
        ch3_pt_file = fullfile(contourDir, ['CH3_AI_pts_3D' suffix '.npy']);
        if(~exist(ch3_pt_file))
            ch3_pt_file = fullfile(contourDir, ['CH3_AI_pts_3D' suffix_ch3 '.npy']);
        end        
        
        if(exist(ch3_pt_file))
            pts = readNPY(ch3_pt_file);
            
            ES_as_first_phase = 1;
            [GLS, MAPSE, V] = compute_GLS(pts, ES_as_first_phase);
            ES_as_first_phase = 0;
            [GLS2, MAPSE2, V] = compute_GLS(pts, ES_as_first_phase);
        
            GLS_3CH = [GLS_3CH; max(GLS2)];
            GLS_3CH_ED_first_phase = [GLS_3CH_ED_first_phase; max(GLS)];
            MAPSE_3CH = [MAPSE_3CH; max(MAPSE2)];
        else
            disp(['-----> cannot find CH3_AI_pts_3D, ' pt_id ' - ' study_date ' - ' suffix]);
            GLS_3CH = [GLS_3CH; -1];
            GLS_3CH_ED_first_phase = [GLS_3CH_ED_first_phase; -1];
            MAPSE_3CH = [MAPSE_3CH; -1];
        end
        
        % --------------------------------------
        
        ch4_pt_file = fullfile(contourDir, ['CH4_AI_pts' suffix '.npy']);
        if(~exist(ch4_pt_file))
            ch4_pt_file = fullfile(contourDir, ['CH4_AI_pts' suffix_ch4 '.npy']);
        end        
        
        if(exist(ch4_pt_file))
            pts = readNPY(ch4_pt_file);
            
            ES_as_first_phase = 1;
            [GLS, MAPSE, lMAPSE] = compute_GLS(pts, ES_as_first_phase);
            ES_as_first_phase = 0;
            [GLS2, MAPSE2, lMAPSE2] = compute_GLS(pts, ES_as_first_phase);
            
            GLS_4CH_2D = [GLS_4CH_2D; max(GLS2)];
            GLS_4CH_ED_first_phase_2D = [GLS_4CH_ED_first_phase_2D; max(GLS)];
            MAPSE_4CH_2D = [MAPSE_4CH_2D; max(MAPSE2)];
            lMAPSE_4CH_2D = [lMAPSE_4CH_2D; max(lMAPSE2)];
        else
            disp(['-----> cannot find CH4_AI_pts, ' pt_id ' - ' study_date ' - ' suffix]);
            GLS_4CH_2D = [GLS_4CH_2D; -1];
            GLS_4CH_ED_first_phase_2D = [GLS_4CH_ED_first_phase_2D; -1];
            MAPSE_4CH_2D = [MAPSE_4CH_2D; -1];
        end
        
        ch2_pt_file = fullfile(contourDir, ['CH2_AI_pts' suffix '.npy']);
        if(~exist(ch2_pt_file))
            ch2_pt_file = fullfile(contourDir, ['CH2_AI_pts' suffix_ch2 '.npy']);
        end        
        
        if(exist(ch2_pt_file))
            pts = readNPY(ch2_pt_file);
            ES_as_first_phase = 1;
            [GLS, MAPSE, V] = compute_GLS(pts, ES_as_first_phase);
            ES_as_first_phase = 0;
            [GLS2, MAPSE2, V] = compute_GLS(pts, ES_as_first_phase);
            
            GLS_2CH_2D = [GLS_2CH_2D; max(GLS2)];
            GLS_2CH_ED_first_phase_2D = [GLS_2CH_ED_first_phase_2D; max(GLS)];
            MAPSE_2CH_2D = [MAPSE_2CH_2D; max(MAPSE2)];
        else
            disp(['-----> cannot find CH2_AI_pts, ' pt_id ' - ' study_date ' - ' suffix]);
            GLS_2CH_2D = [GLS_2CH_2D; -1];
            GLS_2CH_ED_first_phase_2D = [GLS_2CH_ED_first_phase_2D; -1];
            MAPSE_2CH_2D = [MAPSE_2CH_2D; -1];
        end
               
        ch3_pt_file = fullfile(contourDir, ['CH3_AI_pts' suffix '.npy']);
        if(~exist(ch3_pt_file))
            ch3_pt_file = fullfile(contourDir, ['CH3_AI_pts' suffix_ch3 '.npy']);
        end        
        
        if(exist(ch3_pt_file))
            pts = readNPY(ch3_pt_file);
            ES_as_first_phase = 1;
            [GLS, MAPSE, V] = compute_GLS(pts, ES_as_first_phase);
            ES_as_first_phase = 0;
            [GLS2, MAPSE2, V] = compute_GLS(pts, ES_as_first_phase);
        
            GLS_3CH_2D = [GLS_3CH_2D; max(GLS2)];
            GLS_3CH_ED_first_phase_2D = [GLS_3CH_ED_first_phase_2D; max(GLS)];
            MAPSE_3CH_2D = [MAPSE_3CH_2D; max(MAPSE2)];
        else
            disp(['-----> cannot find CH3_AI_pts, ' pt_id ' - ' study_date ' - ' suffix]);
            GLS_3CH_2D = [GLS_3CH_2D; -1];
            GLS_3CH_ED_first_phase_2D = [GLS_3CH_ED_first_phase_2D; -1];
            MAPSE_3CH_2D = [MAPSE_3CH_2D; -1];
        end
    end
    
    record = table(patients, pids, Age, Gender, Height, Weight, BMI, BSA, EF, ESV, EDV, ... 
        MASS, Infarction, NonIschemicScar, LGE, FELKER, ...
        GLS_Strain, ...
        GLS_4CH_ED_first_phase, GLS_3CH_ED_first_phase, GLS_2CH_ED_first_phase, ...
        GLS_4CH, GLS_3CH, GLS_2CH, ...
        GLS_4CH_ED_first_phase_2D, GLS_3CH_ED_first_phase_2D, GLS_2CH_ED_first_phase_2D, ...
        GLS_4CH_2D, GLS_3CH_2D, GLS_2CH_2D, ...
        MAPSE_4CH, MAPSE_3CH, MAPSE_2CH, lMAPSE_4CH, ...
        MAPSE_4CH_2D, MAPSE_3CH_2D, MAPSE_2CH_2D, lMAPSE_4CH_2D, ...
        pt_dirs );
end

function [GLS, MAPSE, lMAPSE] = compute_GLS(pts, ES_as_first_phase)

    if(~isempty(find(pts<2)))
        GLS = -1;
        MAPSE = -1;
        lMAPSE = -1;
        return
    end

    GLS = -1;
    MAPSE = -1;
    
    ptc = 0.5 * (pts(1,:,:) + pts(2,:,:));
    ptc = squeeze(ptc);
    
    N = size(pts, 3);
    
    L = zeros(N,1);
    
    for j = 1:N        
        L(j) = norm([ptc(1,j), ptc(2,j)] - [pts(3,1,j), pts(3,2,j)]);
    end

    if(ES_as_first_phase)
        GLS = 100 * (L(1) - L) / L(1);
    else
        maxL = max(L(:));
        sorted_L = sort(L);
        %GLS = 100 * (maxL - L) / mean(sorted_L(N-5:N));
        % maxL = (0.5 * (L(1)+L(end)));
        GLS = 100 * (maxL - L) / (0.5 * (L(1)+L(end)));
    end
    
    MAPSE = zeros(N, 1);
    ptc_ED = (ptc(:,1) + ptc(:,end)) / 2;
    for j = 1:N        
        MAPSE(j) = norm([ptc(1,j), ptc(2,j)] - [ptc_ED(1,1), ptc_ED(2,1)]);
    end
    
    lMAPSE = zeros(N, 1);
    ptl_ED = squeeze((pts(2, :,1) + pts(2, :, end)) / 2);
    for j = 1:N        
        lMAPSE(j) = norm([pts(2,1,j), pts(2,2,j)] - [ptl_ED(1), ptl_ED(2)]);
    end
end

% function GLS = compute_GLS(pts, ES_as_first_phase)
% 
%     if(~isempty(find(pts<2)))
%         GLS = -1;
%         return
%     end
% 
%     GLS = -1;
%     
%     ptc = 0.5 * (pts(1,:,:) + pts(2,:,:));
%     ptc = squeeze(ptc);
%     
%     N = size(pts, 3);
%     
%     L = zeros(N,1);
%     
%     for j = 1:N        
%         L(j) = norm([ptc(1,j), ptc(2,j)] - [pts(3,1,j), pts(3,2,j)]);
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
% end

