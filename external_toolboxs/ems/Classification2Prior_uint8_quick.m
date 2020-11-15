
function mix = Classification2Prior_uint8_quick(mix)
% transform the classification to prior
% for the 5 class

% prior: csf, cortex, wm1, wm2, outlier
% classification: wm1, wm2, cortex, csf, outlier

ndata = size(mix.indexes, 1);

x = mix.indexes(:, 2);
y = mix.indexes(:, 1);
z = mix.indexes(:, 3);
f = ones(ndata, 1, 'uint32');

ind = sub2ind(size(mix.classification), y, x, z, f(:));
mix.priors(:, 3) = double(mix.classification(ind))/255;

ind = sub2ind(size(mix.classification), y, x, z, 2*f(:));
mix.priors(:, 4) = double(mix.classification(ind))/255;

ind = sub2ind(size(mix.classification), y, x, z, 3*f(:));
mix.priors(:, 2) = double(mix.classification(ind))/255;

ind = sub2ind(size(mix.classification), y, x, z, 4*f(:));
mix.priors(:, 1) = double(mix.classification(ind))/255;

mix.priors(:, 5) = 1 - sum(mix.priors(:, 1:4), 2);
return;