function labels = prepare_SAX_label(data_dir, label_dir, lax_label, plot_flag, max_N)
% labels = prepare_SAX_label(data_dir, lax_label, plot_flag, max_N)
% if no landmarks for a image, landmarks is -1
% SAX can have up to 5 landmarks, RVI-Anterior, RVI-Posterior, LV Center, end-point-myo, start-point_myo
% if only three landmarks are avaible, they are RVI-Anterior, RVI-Posterior, LV Center
% if only one landmark is available, it is LV Center

labels = [];
fields = fieldnames(lax_label)

if(nargin<5)
    max_N = numel(fields);
end

if(max_N>numel(fields))
    max_N = numel(fields);
end

wrong_label_cases = [];

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
    
    disp(['load ' num2str(k) ' out of ' num2str(max_N) ' - ' v_name]);
    
    pt1 = [-1 -1];
    pt2 = [-1 -1];
    pt3 = [-1 -1];
    pt4 = [-1 -1];
    pt5 = [-1 -1];
    
    if(numel(v.regions)==3)
        pt1 = [v.regions(1).shape_attributes.cx v.regions(1).shape_attributes.cy];
        pt2 = [v.regions(2).shape_attributes.cx v.regions(2).shape_attributes.cy];
        pt3 = [v.regions(3).shape_attributes.cx v.regions(3).shape_attributes.cy];
        
        if(start_e1>0)
            pt1(1) = pt1(1) - start_e1 -1;
            pt2(1) = pt2(1) - start_e1 -1;
            pt3(1) = pt3(1) - start_e1 -1;
        end
        
        if(start_ro>0)
            pt1(2) = pt1(2) - start_ro -1;
            pt2(2) = pt2(2) - start_ro -1;
            pt3(2) = pt3(2) - start_ro -1;
        end
        
        if(E1>E1_im)
            r = E1/E1_im;
            pt1(1) = pt1(1)/r;
            pt2(1) = pt2(1)/r;
            pt3(1) = pt3(1)/r;
        end
        
        if(RO>RO_im)
            r = RO/RO_im;
            pt1(2) = pt1(2)/r;
            pt2(2) = pt2(2)/r;
            pt3(2) = pt3(2)/r;
        end
        
        if(plot_flag)
            figure; imagescn(im);
            hold on
            plot(pt1(1)+1, pt1(2)+1, 'ro');
            plot(pt2(1)+1, pt2(2)+1, 'bo');
            plot(pt3(1)+1, pt3(2)+1, 'yo');
            
            pause
            closeall
        end
    elseif(numel(v.regions)==1)
        pt3 = [v.regions(1).shape_attributes.cx v.regions(1).shape_attributes.cy];
        
        if(start_e1>0)
            pt3(1) = pt3(1) - start_e1 -1;
        end
        
        if(start_ro>0)
            pt3(2) = pt3(2) - start_ro -1;
        end
        
        if(E1>E1_im)
            r = E1/E1_im;
            pt3(1) = pt3(1)/r;
        end
        
        if(RO>RO_im)
            r = RO/RO_im;
            pt3(2) = pt3(2)/r;
        end
        
        if(plot_flag)
            figure; imagescn(im);
            hold on
            plot(pt3(1)+1, pt3(2)+1, 'yo');
            
            pause
            closeall
        end
    elseif(numel(v.regions)==5)
        pt1 = [v.regions(1).shape_attributes.cx v.regions(1).shape_attributes.cy];
        pt2 = [v.regions(2).shape_attributes.cx v.regions(2).shape_attributes.cy];
        pt3 = [v.regions(3).shape_attributes.cx v.regions(3).shape_attributes.cy];
        pt4 = [v.regions(4).shape_attributes.cx v.regions(4).shape_attributes.cy];
        pt5 = [v.regions(5).shape_attributes.cx v.regions(5).shape_attributes.cy];
        
        if(start_e1>0)
            pt1(1) = pt1(1) - start_e1 -1;
            pt2(1) = pt2(1) - start_e1 -1;
            pt3(1) = pt3(1) - start_e1 -1;
            pt4(1) = pt4(1) - start_e1 -1;
            pt5(1) = pt5(1) - start_e1 -1;
        end
        
        if(start_ro>0)
            pt1(2) = pt1(2) - start_ro -1;
            pt2(2) = pt2(2) - start_ro -1;
            pt3(2) = pt3(2) - start_ro -1;
            pt4(2) = pt4(2) - start_ro -1;
            pt5(2) = pt5(2) - start_ro -1;
        end
        
        if(E1>E1_im)
            r = E1/E1_im;
            pt1(1) = pt1(1)/r;
            pt2(1) = pt2(1)/r;
            pt3(1) = pt3(1)/r;
            pt4(1) = pt4(1)/r;
            pt5(1) = pt5(1)/r;
        end
        
        if(RO>RO_im)
            r = RO/RO_im;
            pt1(2) = pt1(2)/r;
            pt2(2) = pt2(2)/r;
            pt3(2) = pt3(2)/r;
            pt4(2) = pt4(2)/r;
            pt5(2) = pt5(2)/r;
        end
        
        if(plot_flag)
            figure; imagescn(im);
            hold on
            plot(pt1(1)+1, pt1(2)+1, 'ro');
            plot(pt2(1)+1, pt2(2)+1, 'bo');
            plot(pt3(1)+1, pt3(2)+1, 'yo');
            plot(pt4(1)+1, pt4(2)+1, 'gx');
            plot(pt5(1)+1, pt5(2)+1, 'g+');
            
            pause
            closeall
        end
    elseif(numel(v.regions)==0)
        if(plot_flag)
            figure; imagescn(im);
            pause
            closeall
        end
    else
        warning(['incorrect labels for ' num2str(k) ' - ' v_name ]);
        wrong_label_cases = [wrong_label_cases; {k, v_name}];
        continue;
    end
    
    labels = [labels; {v_name, [pt1;pt2;pt3;pt4;pt5]}];
    
    check_valid_pt(RO_im, E1_im, pt1);
    check_valid_pt(RO_im, E1_im, pt2);
    check_valid_pt(RO_im, E1_im, pt3);
    check_valid_pt(RO_im, E1_im, pt4);
    check_valid_pt(RO_im, E1_im, pt5);
end

wrong_label_cases

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
