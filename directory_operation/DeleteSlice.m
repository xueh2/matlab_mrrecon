
function DeleteSlice(inputName, outputName, slicesNum)
% load a hdr image and delete specific slice
[data, header] = LoadAnalyze(inputName, 'Grey');

pp = sort(slicesNum);

num = numel(pp);
header2 = header;
header2.zsize = header.zsize - num;
dataFinal = zeros([header2.ysize header2.xsize header2.zsize], 'uint32');

index = 0;
dataFinal(:,:,1:pp(1)-1) = data(:,:,1:pp(1)-1);
index = pp(1);
for i=2:num
    dataFinal(:,:,index:index+pp(i)-1-pp(i-1)-1) = data(:,:,pp(i-1)+1:pp(i)-1);
    index = index+pp(i)-1-pp(i-1);
end
dataFinal(:,:,index:header2.zsize) = data(:,:,pp(num)+1:header.zsize);

SaveAnalyze(dataFinal, header2, outputName, 'Grey');