
function mix = Classification2Prior_uint8_4classes_quick(mix)
% transform the classification to prior
% for the 5 class

% prior: csf, cortex, wm, outlier
% classification: wm, cortex, csf, outlier

ndata = size(mix.indexes, 1);

% for i = 1:ndata
%     x = mix.indexes(i, 2);
%     y = mix.indexes(i, 1);
%     z = mix.indexes(i, 3);
%     mix.priors(i, 3) = double(mix.classification(y, x, z, 1))/255;
%     mix.priors(i, 2) = double(mix.classification(y, x, z, 2))/255;
%     mix.priors(i, 1) = double(mix.classification(y, x, z, 3))/255;
%     mix.priors(i, 4) = 1 - sum(mix.priors(i, 1:3));
% end

x = mix.indexes(:, 2);
y = mix.indexes(:, 1);
z = mix.indexes(:, 3);
f = ones(ndata, 1, 'uint32');

ind = sub2ind(size(mix.classification), y, x, z, f(:));
mix.priors(:, 3) = double(mix.classification(ind))/255;

ind = sub2ind(size(mix.classification), y, x, z, 2*f(:));
mix.priors(:, 2) = double(mix.classification(ind))/255;

ind = sub2ind(size(mix.classification), y, x, z, 3*f(:));
mix.priors(:, 1) = double(mix.classification(ind))/255;

mix.priors(:, 4) = 1 - sum(mix.priors(:, 1:3), 2);

return;