function sector_contour = get_sector_contour(sectors, s, new_len, n_components)
% sector_contour = get_sector_contour(sectors, s, new_len, n_components)

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
