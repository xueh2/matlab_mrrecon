
function res = SplitContourForAHAModel(endo, lv_center, rv_insertion, num_seg)
% res = SplitContourForAHAModel(endo, lv_center, rv_insertion, num_seg)
% split contour into num_seg portions

num_pt = size(endo, 1);

%  compute angles for all points
pt_angle = zeros(num_pt, 1);

lv_vec = lv_center - rv_insertion;

res = cell(num_seg, 1);

s_angle = 360/num_seg;

for n=1:num_pt    
    vvec = endo(n, :) - lv_center;
    
    pt_angle(n) = atan2d(lv_vec(1)*vvec(2)-lv_vec(2)*vvec(1), lv_vec(1)*vvec(1)+lv_vec(2)*vvec(2));
    
%     u = [lv_vec(1) lv_vec(2) 0];
%     v = [vvec(1) vvec(2) 0];
%     pt_angle(n) = atan2d(norm(cross(u,v)),dot(u,v));
    if(pt_angle(n)<0)
        pt_angle(n) = abs(pt_angle(n)) + 180;
    end
    
    sector = floor(pt_angle(n) / s_angle) + 1;
    res{sector} = [ res{sector}; endo(n, :)];
end

for s=1:num_seg  
    pt_list = res{s};
    num_pt = size(pt_list, 1);
    if(num_pt==0)
        continue;
    end
    
    pt_angle = zeros(num_pt, 1);
    for n=1:num_pt    
        vvec = pt_list(n, :) - lv_center;

        pt_angle(n) = atan2d(lv_vec(1)*vvec(2)-lv_vec(2)*vvec(1), lv_vec(1)*vvec(1)+lv_vec(2)*vvec(2));   
        if(pt_angle(n)<0)
            pt_angle(n) = abs(pt_angle(n)) + 180;
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

