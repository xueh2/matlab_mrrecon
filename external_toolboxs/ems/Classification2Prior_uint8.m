
function mix = Classification2Prior_uint8(mix)
% transform the classification to prior
% for the 5 class

% prior: csf, cortex, wm1, wm2, outlier
% classification: wm1, wm2, cortex, csf, outlier

ndata = size(mix.indexes, 1);

for i = 1:ndata
    x = mix.indexes(i, 2);
    y = mix.indexes(i, 1);
    z = mix.indexes(i, 3);
    mix.priors(i, 3) = double(mix.classification(y, x, z, 1))/255;
    mix.priors(i, 4) = double(mix.classification(y, x, z, 2))/255;
    mix.priors(i, 2) = double(mix.classification(y, x, z, 3))/255;
    mix.priors(i, 1) = double(mix.classification(y, x, z, 4))/255;
    mix.priors(i, 5) = 1 - sum(mix.priors(i, 1:4));
end

return;