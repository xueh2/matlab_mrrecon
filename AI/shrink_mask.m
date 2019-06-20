
function [endo_mask_res, epi_mask_res] = shrink_mask(endo_mask, epi_mask, ratio_endo, ratio_epi)
% [endo_mask_res, epi_mask_res] = shrink_mask(endo_mask, epi_mask, ratio_endo, ratio_epi)
% if ratio_endo is -0.1, shrink by 10%; if ratio_endo is 0.1, expand by 10%; if ratio_endo is 0, no change

myo_mask = epi_mask - endo_mask;

endo_c = get_sector_contour(endo_mask, 1, 4000);
epi_c = get_sector_contour(epi_mask, 1, 4000);

k_endo = kdtree(endo_c);
[idx,pout] = kdtree_closestpoint(k_endo, epi_c);
wt_epi = epi_c-pout;
wt_epi = sqrt(sum(wt_epi.*wt_epi, 2));

k_epi = kdtree(epi_c);
[idx,pout] = kdtree_closestpoint(k_epi, endo_c);
wt_endo = endo_c-pout;
wt_endo = sqrt(sum(wt_endo.*wt_endo, 2));

[I, J] = find(epi_mask>0);
lv_center = [mean(J) mean(I)];

endo_c_res = shrink_contour(endo_c, ratio_endo, wt_endo, lv_center);
epi_c_res = shrink_contour(epi_c, ratio_epi, wt_epi, lv_center);

figure; imagescn(myo_mask);
hold on
plot(endo_c(:,1), endo_c(:,2), 'y-');
plot(epi_c(:,1), epi_c(:,2), 'g-');
plot(endo_c_res(:,1), endo_c_res(:,2), 'y--');
plot(epi_c_res(:,1), epi_c_res(:,2), 'g--');
hold off

endo_mask_res = roipoly(endo_mask, endo_c_res(:,1), endo_c_res(:,2));
epi_mask_res = roipoly(epi_mask, epi_c_res(:,1), epi_c_res(:,2));

end

function c_res = shrink_contour(c, ratio, wt, lv_center)

c_res = c;
if(ratio==1)
    return;
end

nt = size(c, 1);

c2 = c - repmat(lv_center, [nt 1]);
d = sqrt(sum(c2.*c2, 2));
v2 = c2 ./ repmat(d, [1 2]);
c2 = c2 + repmat(wt, [1 2]) .* (ratio-1) .* v2;
c_res = c2 + repmat(lv_center, [nt 1]);
end
    