
function mix = Prior2Classification_uint8_quick(mix)
% transform the prior to classification
% for the 5 class

% prior: csf, cortex, wm1, wm2, outlier
% classification: wm1, wm2, cortex, csf, outlier

ndata = size(mix.indexes, 1);

x = uint32(mix.indexes(:, 2));
y = uint32(mix.indexes(:, 1));
z = uint32(mix.indexes(:, 3));
f = ones(ndata, 1, 'uint32');

temp = zeros(ndata, 1, 'uint32');

ind = sub2ind(size(mix.classification), y, x, z, f(:));
mix.classification(ind) = uint8(round(255*mix.priors(:, 3)));
temp = temp + uint32(mix.classification(ind));

ind = sub2ind(size(mix.classification), y, x, z, 2*f(:));
mix.classification(ind) = uint8(round(255*mix.priors(:, 4)));
temp = temp + uint32(mix.classification(ind));

ind = sub2ind(size(mix.classification), y, x, z, 3*f(:));
mix.classification(ind) = uint8(round(255*mix.priors(:, 2)));
temp = temp + uint32(mix.classification(ind));

ind = sub2ind(size(mix.classification), y, x, z, 4*f(:));
mix.classification(ind) = uint8(round(255*mix.priors(:, 1)));
temp = temp + uint32(mix.classification(ind));

ind = sub2ind(size(mix.classification), y, x, z, 5*f(:));
mix.classification(ind) = uint8(255 - temp);
return;