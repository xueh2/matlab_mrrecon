
function centres = ClusteringSuspectedPoints(suspected_volume, classnumber, trynumber, initialCentres )
% use the kmeans clustering to get clustering centres

autoIntial = true;
if ( (isempty(initialCentres) == 0) & (size(initialCentres, 1)==classnumber) & (size(initialCentres, 3)==trynumber) )
    autoIntial = false;
end

index = find(suspected_volume==1);
[row, col, depth] = ind2sub(size(suspected_volume), index);
x = [col row depth];% x, y, z

if ( autoIntial )
    [IDX,C,sumd,D] = kmeans(x, classnumber, 'distance', 'sqEuclidean', 'display', 'iter', 'replicates', trynumber, 'Maxiter', 200, 'EmptyAction', 'drop');
else
    [IDX,C,sumd,D] = kmeans(x, classnumber, 'start', initialCentres, 'distance', 'sqEuclidean', 'display', 'iter', 'Maxiter', 200, 'EmptyAction', 'drop');
end

centres = C;

return
