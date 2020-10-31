function labels = prepare_LAX_segmentation_mask(data_dir, label_dir, lax_label, plot_flag, max_N, load_myo_pts)
% labels = prepare_LAX_segmentation_mask(data_dir, lax_label, plot_flag, max_N)
% if no landmarks for a image, landmarks is -1

fields = fieldnames(lax_label)
if(nargin<5)
    max_N = numel(fields);
end

if(nargin<6)
    load_myo_pts = 0;
end

if(max_N>numel(fields))
    max_N=numel(fields);
end

labels = [];
for k=1:max_N
    v = getfield(lax_label, fields{k});
    [v_path, v_name, ext] = fileparts(v.filename);
    try
        im = readNPY(fullfile(data_dir, v_name));
        im_label = imread(fullfile(label_dir, [v_name '.jpg']));
        
        RO_im = size(im, 1);
        E1_im = size(im, 2);
        
        im_label = single(im_label(:,:,1));
        RO = size(im_label, 1);
        E1 = size(im_label, 2);
        
        if(E1>E1_im)
            ind = [];
            for e1=1:E1
                v_pts = im_label(:,e1);
                if(mean(v_pts(:))==255)
                    ind = [ind; e1];
                end
            end
            if(~isempty(ind))
                start_e1 = max(ind(find(ind<E1/2)));
                if(isempty(start_e1))
                    start_e1 = 0;
                end
                end_e1 = min(ind(find(ind>E1/2)));

                im_label = im_label(:, start_e1+1:end_e1-1);
            else
                start_e1 = 0;
            end
        else
            start_e1 = 0;
        end
        
        if(RO>RO_im)
            ind = [];
            for ro=1:RO
                v_pts = im_label(ro,:);
                if(mean(v_pts(:))==255)
                    ind = [ind; ro];
                end
            end
            if(~isempty(ind))
                start_ro = max(ind(find(ind<RO/2)));
                if(isempty(start_ro))
                    start_ro = 0;
                end
                end_ro = min(ind(find(ind>RO/2)));

                im_label = im_label(start_ro+1:end_ro-1, :);
            else
                start_ro = 0;
            end
        else
            start_ro = 0;            
        end
        
        RO = size(im_label, 1);
        E1 = size(im_label, 2);
        
    catch e
        disp(e.message)
        continue
    end
    
    disp(['load ' num2str(k) ' - ' v_name]);
    
    num_regions = numel(v.regions);
    regions = [];
    for k=1:num_regions
        if(strcmp(v.regions(k).shape_attributes.name, 'polyline'))
            pts_x = v.regions(k).shape_attributes.all_points_x;
            pts_y = v.regions(k).shape_attributes.all_points_y;
            
            bw = poly2mask(pts_x,pts_y,size(im,1), size(im,2));
            
            regions = [regions; {k, pts_x, pts_y, bw, im}];
        end
    end
            
    labels = [labels; {v_name, regions}];
end
end
