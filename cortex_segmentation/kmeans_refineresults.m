
function rkmeans_Seg = kmeans_refineresults(kmeans_Seg, R_open, R_close)
% refine the output of kmeans segmentation

% use close process to remove small seperated parts

R = R_open;
se = strel('ball', R, R, 0);
I_opened = imopen(kmeans_Seg,se);

R = R_close;
se = strel('ball', R, R, 0);
rkmeans_Seg = imclose(I_opened, se);

return