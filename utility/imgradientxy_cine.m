function [dx_cine, dy_cine] = imgradientxy_cine(cine, method)

cine_sz = size(cine);
if(numel(cine_sz)==2)
    slices = 1;
else
    slices = cine_sz(3);
end
dx_cine = zeros(cine_sz);
dy_cine = zeros(cine_sz);

for i = 1:slices
    [dx, dy] = imgradientxy(cine(:, :, i), method);
    dx_cine(:, :, i) = dx;
    dy_cine(:, :, i) = dy;
end
