function frame_index = ComputeMoCoPerformance_find_seg_frames(seg2D)
% frame_index: record which frames have the segmented contours

N = size(seg2D, 3);

frame_index = [];
for n=1:N
    s2D = seg2D(:,:,n);
    s2D = double(s2D);
    
    ind = find(s2D(:)>0);
    if (~isempty(ind))
        frame_index = [frame_index; n];
    end
end
