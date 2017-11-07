
function SaveSlice(inputName, outputName, slicesNum)
% load a hdr image and delete specific slice
[data, header] = LoadAnalyze(inputName, 'Grey');

pp = sort(slicesNum);

num = numel(pp);
header2 = header;
header2.zsize = num;
dataFinal = zeros([header2.ysize header2.xsize header2.zsize], 'uint32');

for i=1:num
   dataFinal(:,:,i) = data(:,:, pp(i)); 
end

SaveAnalyze(dataFinal, header2, outputName, 'Grey');