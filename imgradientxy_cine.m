function [dr_cine, de_cine] = imgradientxy_cine(cine, method)

cine_sz = size(cine);
if (numel(cine_sz) == 2)
    slices = 1;
else
    slices = imsz(3);
end
de_cine = zeros(cine_sz);
dr_cine = zeros(cine_sz);

for i = 1:slices
    [dr, de] = imgradientxy(cine(:, :, i), method);
    de_cine(:, :, i) = -de;
    dr_cine(:, :, i) = dr;
end