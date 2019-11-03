
function [sectors, sector_contours, sector_endo_contours, sector_epi_contours] = mask2sectors(endo_mask, epi_mask, rv_mask, rvi_mask, num_sectors)
% sectors is a mask with 1 being the first sector and 2 being the 2nd sector ...

    plot_flag = 0;    
    new_len = 1000;        
    n_components = 50;
    
    RO = size(endo_mask, 1);
    E1 = size(endo_mask, 2);

    % find lv center
    [I, J] = find(endo_mask>0);
    lv_center = [mean(J) mean(I)];
    figure; imagescn(endo_mask); hold on; plot(lv_center(1), lv_center(2), 'b+');

    % find rv center
    [I, J] = find(rv_mask>0);
    rv_center = [mean(J) mean(I)];
    figure; imagescn(rv_mask); hold on; plot(rv_center(1), rv_center(2), 'b+');

    lv_center2 = img_to_xy(lv_center, RO, E1);
    rv_center2 = img_to_xy(rv_center, RO, E1);
    
    % find rvi
    if(isempty(rvi_mask))
        % assume image is in normal orientation
        
        [I, J] = find(rv_mask>0);
        
        num_pt_rv = numel(I);
        
        rv_vec = [rv_center2(1)-lv_center2(1), rv_center2(2)-lv_center2(2)];

        rvi = [];
        max_angle = 0;
        for pt=1:num_pt_rv
            
            rv_pt = img_to_xy([J(pt) I(pt)], RO, E1);
            
            % angle to rv_center - lv_center            
            rv_angle = get_angle(rv_pt-lv_center2,  rv_vec);
            
            if(rv_angle<=180 & rv_angle>max_angle)
                max_angle = rv_angle;
                rvi = [J(pt) I(pt)];
            end
        end
        
        rvi_mask = zeros(size(endo_mask));
        rvi_mask(rvi) = 1;
    else
        [I, J] = find(rvi_mask>0);
        rvi = [mean(J) mean(I)];
        figure; imagescn(rvi_mask); hold on; plot(rvi(1), rvi(2), 'b+');
    end
    
    rvi2 = img_to_xy(rvi, RO, E1);
    
    % split endo/epi to sectors
    rvi_vec = [rvi2(1)-lv_center2(1), rvi2(2)-lv_center2(2)];
    rv_vec = [rv_center2(1)-lv_center2(1), rv_center2(2)-lv_center2(2)];
    rv_rvi_angle = get_angle(rv_vec, rvi_vec);

    delta_rvi_angle = 360/num_sectors;
    disp(['delta_rvi_angle is ' num2str(delta_rvi_angle)]);

    sectors = zeros(RO, E1);

    myo_mask = epi_mask - endo_mask;
    [I, J] = find(myo_mask>0);
    N_myo_pts = size(I, 1);
    angle_myo_pts = zeros(N_myo_pts,1);

    for n = 1:N_myo_pts
        myo_pts_xy = img_to_xy([J(n) I(n)], RO, E1);
        angle_myo_pts(n) = get_angle(rvi_vec, [myo_pts_xy(1)-lv_center2(1), myo_pts_xy(2)-lv_center2(2)]);
        if (rv_rvi_angle>=180) % rotate rvi clock wise 
            angle_myo_pts(n) = 360 - angle_myo_pts(n);
        end
        sector_no = floor(angle_myo_pts(n)/delta_rvi_angle) +1;
        
        if(sector_no==1)
            sectors(I(n), J(n)) = sector_no;
        else
            sectors(I(n), J(n)) = num_sectors+2-sector_no;
        end
    end

    figure;imagescn(sectors+rv_mask+rvi_mask);PerfColorMap;
    
    % get the contours for each sector
    sector_contours = cell(num_sectors,1);
    for n = 1:num_sectors
        sector_contours{n} = get_sector_contour(sectors, n, new_len, n_components);
    end
    
    figure;imagescn(sectors+rv_mask+rvi_mask);
    hold on
    for n = 1:num_sectors
        cc = sector_contours{n};
        if(isempty(cc))
            continue;
        end
        plot(cc(:,1), cc(:,2), '-');
        text(mean(cc(:,1)), mean(cc(:,2)), [num2str(n) ')'], 'FontSize', 18, 'Color', [1 0 0]);
    end
    
    P = mask2poly(logical(endo_mask));
    endo = [P.X' P.Y'];
    
    P = mask2poly(logical(epi_mask));
    try
        epi = [P.X' P.Y'];
    catch
        if(numel(P)>1)
            
            maxLen = P(1).Length;
            ind = 1;
            for k=2:numel(P)
                if(P(k).Length>maxLen)
                    maxLen = P(k).Length;
                    ind = k;
                end
            end
            
            epi = [P(ind).X' P(ind).Y'];
        end
    end
    
    endo_longer = resample_contour(endo, new_len);
    epi_longer = resample_contour(epi, new_len);
    
    fil_endo = filter_contour(endo_longer, n_components);
    fil_epi = filter_contour(epi_longer, n_components);
    
    figure;imagescn(sectors+rv_mask+rvi_mask);
    hold on;
    plot(fil_endo(:,1), fil_endo(:,2), 'r-');
    plot(fil_epi(:,1), fil_epi(:,2), 'g-');    
    plot(lv_center(1), lv_center(2), 'b+');
    plot(rvi(1), rvi(2), 'b+');
    hold off
    
    [res, sector_endo_contours, sector_epi_contours] = SplitEndoEpiContourForAHAModel(fil_endo, fil_epi, lv_center, rv_center, rvi, num_sectors, RO, E1, plot_flag);
    
    figure;imagescn(sectors+rv_mask+rvi_mask);
    hold on    
    try
        for n = 1:num_sectors
            cc = sector_endo_contours{n};
            plot(cc(:,1), cc(:,2), 'b--');
    %         text(mean(cc(:,1)), mean(cc(:,2)), [num2str(n) ')'], 'FontSize', 18, 'Color', [1 0 0]);
            cc = sector_epi_contours{n};
            plot(cc(:,1), cc(:,2), 'r--');
            text(mean(cc(:,1)), mean(cc(:,2)), [num2str(n) ')'], 'FontSize', 18, 'Color', [1 0 0]);
        end
    catch
    end
end

function r = get_angle(a, b)
    % angle from a to b (rotate a to b)
    % positve angle for counter-clock wise
    % 0-360 degrees
    
    v1_theta = atan2(a(2), a(1));
    v2_theta = atan2(b(2), b(1));

    r = (v2_theta - v1_theta) * (180.0 / pi);

    if r < 0
        r = r + 360.0;
    end
end

function pt = img_to_xy(pt, RO, E1)
%     pt = pt -1;
%     a = [pt(2), E1-1-pt(1)];
%     pt = a;

    pt = pt -1;
    a = [pt(1) RO-1-pt(2)];
    pt = a;

end

function sector_contour = get_sector_contour(sectors, s, new_len, n_components)

    if (nargin<3)
        new_len = 0;
    end

    if (nargin<4)
        n_components = 50;
    end
    
    t = zeros(size(sectors));
    ind = find(sectors(:)==s);
    t(ind(:)) = 1;

    P = mask2poly(logical(t)); 
    
    p_ind = 1;
    p_len = 0;
    for pp=1:numel(P)
        if(P(pp).Length>p_len)
            p_ind = pp;
            p_len = P(pp).Length;
        end
    end
    
    if(p_len>0)
    
        sector_contour = [P(p_ind).X' P(p_ind).Y'];

        if(new_len>0)
            sector_contour = resample_contour(sector_contour, new_len);
            sector_contour = filter_contour(sector_contour, n_components);
        end
    else
        sector_contour = [];
    end
%     [c_s, c_e] = CCMS_Contour(t, 0.5, 4, 0);
%     
%     ptN = size(c_s,1);        
%     sector_contour = zeros(2*ptN, 2);
%     for pt=1:ptN
%         sector_contour( 2*(pt-1)+1, :) = c_s(pt, 2:-1:1);
%         sector_contour( 2*(pt-1)+2, :) = c_e(pt, 2:-1:1);
%     end
end

function fil_c = filter_contour(c, n_components)
    n = numel(c(:,1));
    nfilt=n-n_components-1;
    f = fft(c(:,1));
    f(floor(n/2+1-nfilt/2):floor(n/2+nfilt/2)) = 0.0;
    smoothed_contour_x = abs(ifft(f));
    
    f = fft(c(:,2));
    f(floor(n/2+1-nfilt/2):floor(n/2+nfilt/2)) = 0.0;
    smoothed_contour_y = abs(ifft(f));
    
    fil_c = [smoothed_contour_x(:) smoothed_contour_y(:)];
end

function l_c = resample_contour(c, new_len)
    c1x = c(:,1);
    c1y = c(:,2);
    c1x_interp = interp1(linspace(0,1,length(c1x)),c1x, linspace(0,1,new_len),'spline');
    c1y_interp = interp1(linspace(0,1,length(c1y)),c1y, linspace(0,1,new_len),'spline');
    l_c = [c1x_interp(:) c1y_interp(:)];
end