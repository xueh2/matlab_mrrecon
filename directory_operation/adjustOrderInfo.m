function newFileInfo = adjustOrderInfo(fileinfo, index)
newFileInfo = fileinfo;

num = numel(fileinfo);
for k=1:num
    newFileInfo(k) = fileinfo(index(k));
end
