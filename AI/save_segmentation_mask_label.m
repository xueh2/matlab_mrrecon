function save_segmentation_mask_label(label_dir, labels, save_images)
% save_segmentation_mask_label(label_dir, labels, save_images)

N = size(labels, 1);

mkdir(label_dir);

for k=1:N
    v_name = labels{k, 1};
    regions = labels{k, 2};
    for k=1:size(regions, 1)
        bw = regions{k, 4};
        writeNPY(single(bw), fullfile(label_dir, [v_name '_mask_' num2str(k) '.npy']));        
    end
    
    if(save_images)
        im = regions{1, 5};
        writeNPY(single(im), fullfile(label_dir, [v_name '.npy']));
    end
end