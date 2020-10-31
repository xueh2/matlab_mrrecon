
function data = GetDatafromLabel(mix, indexInPrior)
% set all voxels marked by indexInPrior to be 1

xsize = mix.header.xsize;
ysize = mix.header.ysize;
zsize = mix.header.zsize;

data = zeros([ysize, xsize, zsize], 'uint32');

num = length(indexInPrior);

for k = 1:num
    
    i = mix.indexes(indexInPrior(k), 1);
    j = mix.indexes(indexInPrior(k), 2);
    k = mix.indexes(indexInPrior(k), 3);
    data(i, j, k) = 1;
end
return;