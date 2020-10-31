
function [gmlabel, wmlabel] = GetClassLabel_2classes(imagedata, SegResult);

% mean intensities

meanI = zeros(2,1);

for i = 1:2
    
    p = imagedata(find(SegResult==i));
    meanI(i) = mean(p);
    
end

ind = find(meanI==min(meanI));
gmlabel = ind;
meanI(ind) = -1;

ind = find(meanI==max(meanI));
wmlabel = ind;
meanI(ind) = -1;

return;


