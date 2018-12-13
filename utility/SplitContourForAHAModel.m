
function res = SplitContourForAHAModel(endo, lv_center, rv_center, rv_insertion, num_seg, RO, E1)
% res = SplitContourForAHAModel(endo, lv_center, rv_center, rv_insertion, num_seg, RO, E1)
% split contour into num_seg portions

num_pt = size(endo, 1);

%  compute angles for all points
pt_angle = zeros(num_pt, 1);

% lv_vec = lv_center - rv_insertion;
% lv_vec = -lv_vec;
res = cell(num_seg, 1);

lv_center2 = img_to_xy(lv_center, RO, E1);
rv_center2 = img_to_xy(rv_center, RO, E1);
rvi2 = img_to_xy(rv_insertion, RO, E1);

% split endo/epi to sectors
rvi_vec = [rvi2(1)-lv_center2(1), rvi2(2)-lv_center2(2)];
rv_vec = [rv_center2(1)-lv_center2(1), rv_center2(2)-lv_center2(2)];
rv_rvi_angle = get_angle(rv_vec, rvi_vec);

s_angle = 360/num_seg;
for n=1:num_pt    
%     vvec = endo(n, :) - lv_center;
%     
%     pt_angle(n) = -atan2d(lv_vec(1)*vvec(2)-lv_vec(2)*vvec(1), lv_vec(1)*vvec(1)+lv_vec(2)*vvec(2));
%     
% %     u = [lv_vec(1) lv_vec(2) 0];
% %     v = [vvec(1) vvec(2) 0];
% %     pt_angle(n) = atan2d(norm(cross(u,v)),dot(u,v));
%     if(pt_angle(n)<0)
% %         pt_angle(n) = abs(pt_angle(n)) + 180;
%         pt_angle(n) = pt_angle(n) + 360;
%     end
% 
%     sector(n) = floor(pt_angle(n) / s_angle) + 1;
%     sector(n) = mod(sector(n),num_seg)+1; % sector 2 (0 s_angle)....
%     res{sector(n)} = [ res{sector(n)}; endo(n, :)];

    myo_pts_xy = img_to_xy([endo(n,1) endo(n,2)], RO, E1);
     pt_angle(n) = get_angle(rvi_vec, [myo_pts_xy(1)-lv_center2(1), myo_pts_xy(2)-lv_center2(2)]);
    if (rv_rvi_angle>=180) % rotate rvi clock wise 
         pt_angle(n) = 360 -  pt_angle(n);
    end
    sector_no = floor( pt_angle(n)/s_angle) +1;

    if(sector_no==1)
        sector(n) = sector_no;
    else
        sector(n) = num_seg+2-sector_no;
    end
    res{sector(n)} = [ res{sector(n)}; endo(n, :)];
end
   
for s=1:num_seg  
    pt_list = res{s};
    num_pt = size(pt_list, 1);
    if(num_pt==0)
        continue;
    end
    
    pt_angle = zeros(num_pt, 1);
    for n=1:num_pt    
%         vvec = pt_list(n, :) - lv_center;
% 
%         pt_angle(n) = -atan2d(lv_vec(1)*vvec(2)-lv_vec(2)*vvec(1), lv_vec(1)*vvec(1)+lv_vec(2)*vvec(2));   
%         if(pt_angle(n)<0)
%             pt_angle(n) = pt_angle(n) + 360;
%         end
    
        myo_pts_xy = img_to_xy([pt_list(n, 1) pt_list(n, 2)], RO, E1);
        pt_angle(n) = get_angle(rvi_vec, [myo_pts_xy(1)-lv_center2(1), myo_pts_xy(2)-lv_center2(2)]);
        if (rv_rvi_angle>=180) % rotate rvi clock wise 
             pt_angle(n) = 360 -  pt_angle(n);
        end
%         u = [lv_vec(1) lv_vec(2) 0];
%         v = [vvec(1) vvec(2) 0];
%         pt_angle(n) = atan2d(norm(cross(u,v)),dot(u,v));
%         if(pt_angle(n)<0)
%             pt_angle(n) = pt_angle(n) + 180;
%         end
    end

    [pt_angle_sorted, inds] = sort(pt_angle);
    
    res{s} = pt_list(inds(:), :);
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
    pt = pt -1;
    a = [pt(2), E1-1-pt(1)];
    pt = a;
end

