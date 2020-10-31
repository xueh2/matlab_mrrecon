
function result = minDist_NOSelction(dists, pixelThres, voxelsize)
% mean, std, max, min, mean(abs), std(abs), max(abs), min(abs), percentage larger than pixelThres

num = length(pixelThres);

N = length(dists);

result = zeros(1, 4+num);

result(1) = mean(dists);
result(2) = std(dists);
result(3) = max(dists);
result(4) = min(dists);

result(5) = mean(abs(dists));
result(6) = std(abs(dists));
result(7) = max(abs(dists));
result(8) = min(abs(dists));

for i = 1:num
    
    distThres = pixelThres(i) * voxelsize;
    
    numBigger = length(find(abs(dists)>distThres));
    percentage = numBigger / N;
        
    result(8+i) = percentage;
end
return;
