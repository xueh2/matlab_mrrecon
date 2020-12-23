function prepare_ai_perf_pd(dataDir, resDir, aiDir, pt_ids, files_record_picked, checkprocessed)
% prepare_ai_perf_pd(dataDir, resDir, aiDir, pt_ids, files_record_picked, checkprocessed)

for pt=1:numel(pt_ids)
    
    closeall
    
    pt_id = pt_ids{pt};
    
    if(isnumeric(pt_id))
        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' num2str(pt_id)]);
    else
        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);
    end
    
    case_stress = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'stress', 1);
    case_rest = PerformGadgetronRecon_SavedIsmrmrd_search_patient(files_record_picked, pt_id, 'rest', 1);
       
    if(numel(case_stress)==0 & numel(case_rest)==0)
       continue;
    end
    if(numel(case_stress)>0)
        study_date = case_stress.study_dates(end,:);
        pt_id = case_stress.patientIDs{end,:};
    end
    if(numel(case_rest)>0)
        study_date = case_rest.study_dates(end,:);
        pt_id = case_rest.patientIDs{end,:};
    end
    disp([pt_id ' - ' study_date]);
            
    dst_dir = fullfile(aiDir, study_date, pt_id);
    if(~exist(dst_dir))
        dst_dir = fullfile(aiDir, pt_id);
    end
    
    % stress
    dst_dir_stress = fullfile(dst_dir, 'stress');
               
    % rest
    dst_dir_rest = fullfile(dst_dir, 'rest');
    
    if(~checkprocessed | ~exist(fullfile(dst_dir_stress, 'pd_for_seg.jpg')))
        try
            [pd_s, pd_moco_s] = process_stress_rest(dataDir, resDir, study_date, case_stress);

            writeNPY(single(pd_s(:,:,1,:)), fullfile(dst_dir_stress, 'pd_for_seg.npy'));
            writeNPY(single(pd_moco_s(:,:,1,:)), fullfile(dst_dir_stress, 'pd_moco_for_seg.npy'));

            h = figure; 
            imagescn(pd_s(:,:,1,:), [], [], [12]);
            saveas(h, fullfile(dst_dir_stress, 'pd_for_seg.jpg'), 'jpg');
        catch
            disp(['--> failed for ' case_stress.file_names{1}]);
        end
    end
    
    if(~checkprocessed | ~exist(fullfile(dst_dir_rest, 'pd_for_seg.jpg')))
        try
            [pd_r, pd_moco_r] = process_stress_rest(dataDir, resDir, study_date, case_rest);

            writeNPY(single(pd_r(:,:,1,:)), fullfile(dst_dir_rest, 'pd_for_seg.npy'));
            writeNPY(single(pd_moco_r(:,:,1,:)), fullfile(dst_dir_rest, 'pd_moco_for_seg.npy'));

            h = figure; 
            imagescn(pd_r(:,:,1,:), [], [], [12]);
            saveas(h, fullfile(dst_dir_rest, 'pd_for_seg.jpg'), 'jpg');
        catch
            disp(['--> failed for ' case_rest.file_names{1}]);
        end
    end
    
    closeall
end

end

function [pd, pd_moco] = process_stress_rest(dataDir, resDir, study_date, case_stress)

    stress_dir = fullfile(resDir, study_date, case_stress.file_names{1});

    h5_file = fullfile(dataDir, study_date, [case_stress.file_names{1} '.h5']);
    dset = ismrmrd.Dataset(h5_file, 'dataset');
    hdr = ismrmrd.xml.deserialize(dset.readxml);
    SLC = hdr.encoding(1).encodingLimits.slice.maximum+1;

    pd_0 = analyze75read(fullfile(stress_dir, 'DebugOutput', 'PD_for_moco_row0'));
    RO = size(pd_0, 1);
    E1 = size(pd_0, 2);

    pd = zeros(RO, E1, 4, SLC);
    dx_pd = zeros(RO, E1, 4, SLC);
    dy_pd = zeros(RO, E1, 4, SLC);

    pd_moco = pd;

    for slc=0:SLC-1
        im = analyze75read(fullfile(stress_dir, 'DebugOutput', ['PD_for_moco_row' num2str(slc)]));
        dx = analyze75read(fullfile(stress_dir, 'DebugOutput', ['deformation_field_x_' num2str(slc)]));
        dy = analyze75read(fullfile(stress_dir, 'DebugOutput', ['deformation_field_y_' num2str(slc)]));
        if(size(im,1)~=RO)
            im = permute(im, [2, 1, 3]);
            dx = permute(dx, [2, 1, 3]);
            dy = permute(dy, [2, 1, 3]);
        end
        pd(:,:,:,slc+1) = im;
        dx_pd(:,:,:,slc+1) = dx(:,:,1:4);
        dy_pd(:,:,:,slc+1) = dy(:,:,1:4);

        warpped = Matlab_gt_apply_deformation_field_reg_2D_series(double(im), dx, dy);
        pd_moco(:,:,:,slc+1) = warpped;
    end

    pd = permute(pd, [2 1 3 4]);
    pd_moco = permute(pd_moco, [2 1 3 4]);    
end
