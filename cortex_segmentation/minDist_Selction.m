
function result = minDist_Selction(dists, dists_SDF, pixelThres, voxelsize)

num = length(pixelThres);

N = length(dists);

result = zeros(1, 4+num);

result(1) = min( [mean(dists) mean(dists_SDF)] );
result(2) = min( [std(dists) std(dists_SDF)] );
result(3) = min( [mean(abs(dists)) mean(abs(dists_SDF))] );
result(4) = min( [std(abs(dists)) std(abs(dists_SDF))] );

for i = 1:num
    
    distThres = pixelThres(i) * voxelsize;
    
    numBigger = length(find(abs(dists)>distThres));
    percentage = numBigger / N;
    
    numBigger = length(find(abs(dists_SDF)>distThres));
    percentage_SDF = numBigger / N;
    
    result(4+i) = min(percentage, percentage_SDF);
end
return;
