function L = make_GLS_movie(im, pts, movie_file)
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