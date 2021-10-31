function labels = prepare_CH3_MV_label(data_dir, label_dir, lax_label, plot_flag, max_N)
% labels = prepare_CH3_MV_label(data_dir, label_dir, lax_label, plot_flag, max_N)
% if no landmarks for a image, landmarks is -1

fields = fieldnames(lax_label)
if(nargin<5)
    max_N = numel(fields);
end

if(max_N>numel(fields))
    max_N=numel(fields);
end

checked_cases = [];

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
    
    pt1 = [-1 -1];
    
    N_pt = numel(v.regions);
    
    if(N_pt<1)
        disp(['---> check ' num2str(k) ' - ' v_name]);
        checked_cases = [checked_cases; {k, v_name}];
        continue;
    end
    
    pt1 = [v.regions(end).shape_attributes.cx v.regions(end).shape_attributes.cy];

    if(start_e1>0)
        pt1(1) = pt1(1) - start_e1 -1;
    end

    if(start_ro>0)
        pt1(2) = pt1(2) - start_ro -1;
    end

    if(E1>E1_im)
        r = E1/E1_im;
        pt1(1) = pt1(1)/r;
    end

    if(RO>RO_im)
        r = RO/RO_im;
        pt1(2) = pt1(2)/r;
    end

    if(plot_flag)
        figure; imagescn(im);
        hold on
        plot(pt1(1)+1, pt1(2)+1, 'ro', 'MarkerSize', 12);

        pause
        closeall
    end
    
    labels = [labels; {v_name, [pt1]}];
    
    check_valid_pt(RO_im, E1_im, pt1);
end

checked_cases

end

function check_valid_pt(RO, E1, pt)
    if(pt(2)>1)
        if(pt(2)>RO)
            error('incorrect pt(2)');
        end
    end
    if(pt(1)>1)
        if(pt(1)>E1)
            error('incorrect pt(1)');
        end
    end
end