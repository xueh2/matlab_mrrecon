function pts = generateTestPoints(spatial, temporal)
% generate testing points from spatial and tempoal locations

t = temporal;
s = spatial;

N = numel(t)*size(s, 1);
pts = zeros(N, 3);

ind = 1;
for jj=1:size(s, 1)
    for ii=1:numel(t)
        pts(ind, :) = [s(jj,:) t(ii)];
        ind = ind + 1;
    end
end