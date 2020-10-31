function [L, L_3D] = make_GLS_video_view(im_file, pts_file, movie_file)
    L = [];
    L_3D = [];
    pts = readNPY(pts_file);
    im = readNPY(im_file);
    L = make_GLS_movie(im, pts, movie_file);
end