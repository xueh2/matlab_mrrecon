
function mix = Classification2Prior_uint8_2classes_quick(mix)
% transform the classification to prior
% for the 3 class

% prior: cortex, wm1, wm2
% classification: wm1, wm2, cortex

ndata = size(mix.indexes, 1);

x = mix.indexes(:, 2);
y = mix.indexes(:, 1);
z = mix.indexes(:, 3);
f = ones(ndata, 1, 'uint32');

ind = sub2ind(size(mix.classification), y, x, z, f(:));
mix.priors(:, 2) = double(mix.classification(ind))/255;

ind = sub2ind(size(mix.classification), y, x, z, 2*f(:));
mix.priors(:, 1) = 1 - sum(mix.priors(:, 2), 2);
return;