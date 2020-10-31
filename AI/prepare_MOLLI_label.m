function [ch4_labels, ch2_labels, ch3_labels, sax_labels] = prepare_MOLLI_label(data_dir, label_dir, label, plot_flag, min_N, max_N, load_myo_pts)
% [ch4_labels, ch2_labels, ch3_labels, sax_labels] = prepare_MOLLI_label(data_dir, label_dir, lax_label, plot_flag, max_N, load_myo_pts)
% if no landmarks for a image, landmarks is -1

fields = fieldnames(label)

if(nargin<5)
    min_N = 1;
end

if(nargin<6)
    max_N = numel(fields);
end

if(nargin<7)
    load_myo_pts = 0;
end

ch4_labels = [];
ch2_labels = [];
ch3_labels = [];
sax_labels = [];

if(max_N>numel(fields))
    max_N = numel(fields);
end

for k=min_N:max_N
    v = getfield(label, fields{k});
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
    
    CMR_view = v.file_attributes.CMRView;

    disp(['load - ' CMR_view ' - ' v_name]);

    if(CMR_view=='CH4')
        ch4_labels = fill_lax_labels(v_name, ch4_labels, v, load_myo_pts, im, im_label, start_ro, start_e1, plot_flag);
    end
    
    if(CMR_view=='CH2')
        ch2_labels = fill_lax_labels(v_name, ch2_labels, v, load_myo_pts, im, im_label, start_ro, start_e1, plot_flag);
    end
    
    if(CMR_view=='CH3')
        ch3_labels = fill_lax_labels(v_name, ch3_labels, v, load_myo_pts, im, im_label, start_ro, start_e1, plot_flag);
    end
    
    if(CMR_view=='SAX')
        sax_labels = fill_sax_labels(v_name, sax_labels, v, im, im_label, start_ro, start_e1, plot_flag);
    end
end

end

function labels = fill_lax_labels(v_name, labels, v, load_myo_pts, im, im_label, start_ro, start_e1, plot_flag)

    RO = size(im_label, 1);
    E1 = size(im_label, 2);
    
    RO_im = size(im, 1);
    E1_im = size(im, 2);

    pt1 = [-1 -1];
    pt2 = [-1 -1];
    pt3 = [-1 -1];
    
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
    end
    
    if(numel(v.regions)==4)
        if(load_myo_pts)
            pt1 = [v.regions(4).shape_attributes.cx v.regions(4).shape_attributes.cy];
        else
            pt1 = [v.regions(1).shape_attributes.cx v.regions(1).shape_attributes.cy];
        end
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
            v = E1/E1_im;
            pt1(1) = pt1(1)/v;
            pt2(1) = pt2(1)/v;
            pt3(1) = pt3(1)/v;
        end
        
        if(RO>RO_im)
            v = RO/RO_im;
            pt1(2) = pt1(2)/v;
            pt2(2) = pt2(2)/v;
            pt3(2) = pt3(2)/v;
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
    end
    
    labels = [labels; {v_name, [pt1;pt2; pt3]}];
end

function labels = fill_sax_labels(v_name, labels, v, im, im_label, start_ro, start_e1, plot_flag)

    RO = size(im_label, 1);
    E1 = size(im_label, 2);
    
    RO_im = size(im, 1);
    E1_im = size(im, 2);

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
        warning(['incorrect labels for ' v_name ]);
        figure; imagescn(im);
        return
    end
    
    labels = [labels; {v_name, [pt1;pt2;pt3;pt4;pt5]}];
end