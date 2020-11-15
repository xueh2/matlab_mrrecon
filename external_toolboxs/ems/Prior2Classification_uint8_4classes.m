
function mix = Prior2Classification_uint8_4classes(mix)
% transform the prior to classification
% for the 5 class

% prior: csf, cortex, wm, outlier
% classification: wm, cortex, csf, outlier

ndata = size(mix.indexes, 1);

% tt = 0;
for i = 1:ndata
    x = mix.indexes(i, 2);
    y = mix.indexes(i, 1);
    z = mix.indexes(i, 3);
    mix.classification(y, x, z, 1) = uint8(round(255*mix.priors(i, 3)));
    mix.classification(y, x, z, 2) = uint8(round(255*mix.priors(i, 2)));
    mix.classification(y, x, z, 3) = uint8(round(255*mix.priors(i, 1)));
    mix.classification(y, x, z, 4) = uint8(255 - sum(mix.classification(y, x, z, 1:3)));
end

return;