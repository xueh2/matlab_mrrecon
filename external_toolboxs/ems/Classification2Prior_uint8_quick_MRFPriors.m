
function priors = Classification2Prior_uint8_quick_MRFPriors(classification, mix, priors)
% transform the classification to prior
% for the 5 class

% prior: csf, cortex, wm1, wm2, outlier
% classification: wm1, wm2, cortex, csf, outlier

ndata = size(mix.indexes, 1);

x = mix.indexes(:, 2);
y = mix.indexes(:, 1);
z = mix.indexes(:, 3);
f = ones(ndata, 1, 'uint32');

ind = sub2ind(size(classification), y, x, z, f(:));
priors(:, 3) = double(classification(ind))/255;

ind = sub2ind(size(classification), y, x, z, 2*f(:));
priors(:, 4) = double(classification(ind))/255;

ind = sub2ind(size(classification), y, x, z, 3*f(:));
priors(:, 2) = double(classification(ind))/255;

ind = sub2ind(size(classification), y, x, z, 4*f(:));
priors(:, 1) = double(classification(ind))/255;

priors(:, 5) = 1 - sum(priors(:, 1:4), 2);
return;