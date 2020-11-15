
function mix = Prior2Classification_uint8_2classes_quick(mix)
% transform the prior to classification
% for the 3 class

% prior: cortex, wm
% classification: wm, cortex

ndata = size(mix.indexes, 1);

x = uint32(mix.indexes(:, 2));
y = uint32(mix.indexes(:, 1));
z = uint32(mix.indexes(:, 3));
f = ones(ndata, 1, 'uint32');

temp = zeros(ndata, 1, 'uint32');

ind = sub2ind(size(mix.classification), y, x, z, f(:));
mix.classification(ind) = uint8(round(255*mix.priors(:, 2)));
temp = temp + uint32(mix.classification(ind));

ind = sub2ind(size(mix.classification), y, x, z, 2*f(:));
mix.classification(ind) = uint8(255 - temp);

return;