
function [gmlabel, wmlabel1, wmlabel2] = GetClassLabel_3classes(imagedata, SegResult);

% mean intensities

meanI = zeros(3,1);

for i = 1:3
    
    p = imagedata(find(SegResult==i));
    meanI(i) = mean(p);
    
end

ind = find(meanI==min(meanI));
gmlabel = ind;
meanI(ind) = -1;

tempI = meanI(find(meanI~=-1));
ind = find(meanI==min(tempI));
wmlabel1 = ind;
meanI(ind) = -1;

ind = find(meanI==max(meanI));
wmlabel2 = ind;
meanI(ind) = -1;

return;


