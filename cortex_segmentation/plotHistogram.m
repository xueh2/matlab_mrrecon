
function plotHistogram(x)
% plot histogram of data x

index = find(x>0);
p = x(index(:));
maxP = max(p);

figure;
hist(p, maxP);
return