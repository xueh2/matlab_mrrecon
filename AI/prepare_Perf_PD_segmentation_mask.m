function labels = prepare_Perf_PD_segmentation_mask(data_dir, label_dir, RAW, plot_flag, max_N)
% labels = prepare_Perf_PD_segmentation_mask(data_dir, label_dir, RAW, plot_flag, max_N)
% RAW is the string read from csv record

RAW_N = size(RAW,1)-1;

if(nargin<5)
    max_N = RAW_N;
end

if(max_N>RAW_N)
    max_N=RAW_N;
end

labels = [];
k=2;

while(k<=max_N)
    
    fname = RAW{k, 1}    
    [v_path, v_name, ext] = fileparts(fname);
    
    ind = strfind(v_name, '_pd_1st');
    f_name = v_name(1:ind-1);
    slc = str2num(v_name(end));
    
    disp([' --> process ' f_name ' - slc '  num2str(slc)]);
    
    try
        im = readNPY(fullfile(data_dir, f_name, 'pd.npy'));
        im_label = imread(fullfile(label_dir, [v_name '.jpg']));
        
        im = im(:,:,1, slc);
        
        % note a permute
        im = im';
        
        im_label = rgb2gray(im_label);
        
        RO_im = size(im, 1);
        E1_im = size(im, 2);
        
        im_label = single(im_label);
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
        
        ratio_RO = RO/RO_im;
        ratio_E1 = E1/E1_im;
        
    catch e
        disp(e.message)
        continue
    end
    
    disp(['load ' num2str(k) ' - ' v_name]);
    
    num_regions = RAW{k,4};
    if(num_regions>0)
        
        regions = [];
        for rr=1:num_regions
            region_label = jsondecode(RAW{k+rr-1, 6});
            
            if(strcmp(region_label.name, 'polyline') | strcmp(region_label.name, 'polygon'))
                pts_x = (region_label.all_points_x-start_e1) / ratio_E1;
                pts_y = (region_label.all_points_y-start_ro) / ratio_RO;

                bw = poly2mask(pts_x,pts_y,size(im,1), size(im,2));

                regions = [regions; {k, pts_x, pts_y, bw, im}];
            end
        end

        labels = [labels; {v_name, regions}];
        
        k = k+num_regions;
        
        if(plot_flag)
           newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
         
            h = figure; imagescn(im, [0 4*mean(im(:))], [], [], 12);
            hold on
            for rr=1:num_regions
                C = mask2contour(regions{rr, 4}, 1, 2000, 24);
                color_ind = mode(rr, 4);
                if(rr==4)
                    color_ind = 1;
                end
                plot(C(:,1), C(:,2), 'Color', newcolors(color_ind, :));
            end
            hold off
            
            pause
            closeall
        end
    else
        k = k+1;
    end
end
end
