function plot_results_ai_landmark_detection_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, do_2D, do_3D, check_processed, suffix)
% plot_results_ai_landmark_detection_retro_binning_rt(resDir, aiDir, pt_ids, files_record_picked, case_4chs, case_2chs, case_3chs, case_saxs, do_2D, do_3D, check_processed, suffix)

    for pt=1:numel(pt_ids)

        closeall

        pt_id = pt_ids{pt};

        disp([num2str(pt) ' out of ' num2str(numel(pt_ids)) ' - ' pt_id]);

        if(numel(pt_ids)==numel(case_4chs))
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

        for d=1:size(case_4ch,1)
            disp(['--> ' case_4ch.file_names{d}]);
            case_4ch_dir = fullfile(resDir, case_4ch.file_names{d});
            [path, sname, ext] = fileparts(case_4ch_dir);
            sname= sname(~isspace(sname));
            dst_dir_4ch = fullfile(dst_dir, ['ch4_' sname]);
            [L, L_3D] = make_video_view(do_2D, do_3D, 'ch4', ['_' sname '_' suffix], check_processed, contourDir, dst_dir, dst_dir_4ch);
        end
        
        for d=1:size(case_3ch,1)
            disp(['--> ' case_3ch.file_names{d}]);
            case_3ch_dir = fullfile(resDir, case_3ch.file_names{d});
            [path, sname, ext] = fileparts(case_3ch_dir);
            sname= sname(~isspace(sname));
            dst_dir_3ch = fullfile(dst_dir, ['ch3_' sname]);
            [L, L_3D] = make_video_view(do_2D, do_3D, 'ch3', ['_' sname '_' suffix], check_processed, contourDir, dst_dir, dst_dir_3ch);
        end
        
        for d=1:size(case_2ch,1)
            disp(['--> ' case_2ch.file_names{d}]);
            case_2ch_dir = fullfile(resDir, case_2ch.file_names{d});
            [path, sname, ext] = fileparts(case_2ch_dir);
            sname= sname(~isspace(sname));
            dst_dir_2ch = fullfile(dst_dir, ['ch2_' sname]);
            [L, L_3D] = make_video_view(do_2D, do_3D, 'ch2', ['_' sname '_' suffix], check_processed, contourDir, dst_dir, dst_dir_2ch);
        end
        
%         if(exist(fullfile(contourDir, 'SAX_AI_probs.npy')))
%             dd = fullfile(dst_dir, 'sax');
%             probs = readNPY(fullfile(contourDir, 'SAX_AI_probs.npy'));
%             pts = readNPY(fullfile(contourDir, 'SAX_AI_pts.npy'));
%             im = readNPY(fullfile(dst_dir_sax, 'data.npy'));
%         end    
    end
end

function L = make_movie(im, pts, movie_file)
    N = size(im, 3);
    if(N<10)
        L = 0;
        return;
    end
    vidObj = VideoWriter(movie_file);
    open(vidObj);
    
    im2 = Matlab_gt_resize_2D_image(double(im), 4*size(im,1), 4*size(im,2), 5);
    pts2 = 4 * pts + 1;
    
    ptc = 0.5 * (pts2(1,:,:) + pts2(2,:,:));
    ptc = squeeze(ptc);
    
    % dataScaled = normalizeWindowSetting(im2, 1.25*median(im2(:)), 3*median(im2(:)));
    dataScaled = normalizeWindowSetting(im2, 3*mean(im2(:)), 4*mean(im2(:)));
    
    L = zeros(N,1);
    
    figure;
    imshow(double(dataScaled(:,:,1)/255));
    hold on
    ax = gca;    
    for j = 1:N
        axes(ax)
        imshow(double(dataScaled(:,:,j)/255), 'Parent', ax);
        plot(pts2(:,1,j), pts2(:,2,j), 'r+', 'MarkerSize', 24);
        plot(ptc(1,j), ptc(2,j), 'b.', 'MarkerSize', 18);
        plot([ptc(1,j), pts2(3,1,j)], [ptc(2,j), pts2(3,2,j)], 'y--', 'LineWidth', 2.0);
        drawnow;
        currFrame = getframe(gca);
        writeVideo(vidObj,currFrame);
        
        L(j) = norm([ptc(1,j), ptc(2,j)] - [pts2(3,1,j), pts2(3,2,j)]);
    end
    close(vidObj);
    
    [path, name, ext] = fileparts(movie_file);
    h = figure; 
    plot([1:N], L);
    xlabel('Phase');
    ylabel('LAX length, mm');
    box on
    saveas(h, fullfile(path, [name '_L.fig']), 'fig');
    saveas(h, fullfile(path, [name '_L.jpg']), 'jpg');
    closeall
end

function [L, L_3D] = make_video_view(do_2D, do_3D, view_str, suffix, check_processed, contourDir, dst_dir, dst_dir_view)

    view_str_upper = upper(view_str);

    L = [];
    L_3D = [];
    
    if(do_2D & exist(fullfile(contourDir, [view_str_upper '_AI_probs' suffix '.npy'])))
        if(~check_processed | ~exist(fullfile(dst_dir_view, ['LV_length' suffix '.npy'])))
            dd = fullfile(dst_dir, view_str);
            probs = readNPY(fullfile(contourDir, [view_str_upper '_AI_probs' suffix '.npy']));
            pts = readNPY(fullfile(contourDir, [view_str_upper '_AI_pts' suffix '.npy']));
            im = readNPY(fullfile(dst_dir_view, 'data.npy'));
            L = make_movie(im, pts, fullfile(dst_dir_view, [view_str_upper '_pts' suffix '.avi']));
            if(numel(L)>10)
                copyfile(fullfile(dst_dir_view, [view_str_upper '_pts' suffix '.avi']), contourDir);
                copyfile(fullfile(dst_dir_view, [view_str_upper '_pts' suffix '_L.fig']), contourDir);
                copyfile(fullfile(dst_dir_view, [view_str_upper '_pts' suffix '_L.jpg']), contourDir);
                writeNPY(single(L), fullfile(dst_dir_view, ['LV_length' suffix '.npy']));
            end
        end
    end        
    if(do_3D & exist(fullfile(contourDir, [view_str_upper '_AI_probs_3D' suffix '.npy'])))
        if(~check_processed | ~exist(fullfile(dst_dir_view, ['LV_length_3D' suffix '.npy'])))
            dd = fullfile(dst_dir, view_str);
            probs = readNPY(fullfile(contourDir, [view_str_upper '_AI_probs_3D' suffix '.npy']));
            pts = readNPY(fullfile(contourDir, [view_str_upper '_AI_pts_3D' suffix '.npy']));
            im = readNPY(fullfile(dst_dir_view, 'data.npy'));
            L = make_movie(im, pts, fullfile(dst_dir_view, [view_str_upper '_pts_3D' suffix '.avi']));
            if(numel(L)>10)
                copyfile(fullfile(dst_dir_view, [view_str_upper '_pts_3D' suffix '.avi']), contourDir);
                copyfile(fullfile(dst_dir_view, [view_str_upper '_pts_3D' suffix '_L.fig']), contourDir);
                copyfile(fullfile(dst_dir_view, [view_str_upper '_pts_3D' suffix '_L.jpg']), contourDir);
                writeNPY(single(L), fullfile(dst_dir_view, ['LV_length_3D' suffix '.npy']));
            end
        end
    end
end