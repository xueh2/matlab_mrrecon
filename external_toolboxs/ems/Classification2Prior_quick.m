
function mix = Classification2Prior_quick(mix)
% transform the classification to prior
% for the 5 class

% prior: csf, cortex, wm1, wm2, outlier
% classification: wm1, wm2, cortex, csf, outlier

ndata = size(mix.indexes, 1);

% for i = 1:ndata
%     mix.priors(i, 3) = mix.classification(mix.indexes(i, 1), mix.indexes(i, 2), mix.indexes(i, 3), 1);
%     mix.priors(i, 4) = mix.classification(mix.indexes(i, 1), mix.indexes(i, 2), mix.indexes(i, 3), 2);
%     mix.priors(i, 2) = mix.classification(mix.indexes(i, 1), mix.indexes(i, 2), mix.indexes(i, 3), 3);
%     mix.priors(i, 1) = mix.classification(mix.indexes(i, 1), mix.indexes(i, 2), mix.indexes(i, 3), 4);
%     mix.priors(i, 5) = mix.classification(mix.indexes(i, 1), mix.indexes(i, 2), mix.indexes(i, 3), 5);
% end

% mix.priors(:, 3) = mix.classification(mix.sampleInd * 1);
% mix.priors(:, 4) = mix.classification(mix.sampleInd * 2);
% mix.priors(:, 2) = mix.classification(mix.sampleInd * 3);
% mix.priors(:, 1) = mix.classification(mix.sampleInd * 4);
% mix.priors(:, 5) = mix.classification(mix.sampleInd * 5);

x = uint32(mix.indexes(:, 2));
y = uint32(mix.indexes(:, 1));
z = uint32(mix.indexes(:, 3));
f = ones(ndata, 1, 'uint32');

ind = sub2ind(size(mix.classification), y, x, z, f(:));
mix.priors(:, 3) = mix.classification(ind);

ind = sub2ind(size(mix.classification), y, x, z, 2*f(:));
mix.priors(:, 4) = mix.classification(ind);

ind = sub2ind(size(mix.classification), y, x, z, 3*f(:));
mix.priors(:, 2) = mix.classification(ind);

ind = sub2ind(size(mix.classification), y, x, z, 4*f(:));
mix.priors(:, 1) = mix.classification(ind);

mix.priors(:, 5) = 1 - sum(mix.priors(:, 1:4), 2);
return;