
function label3D = GetEachRegion(labelvolume, label)


label3D = zeros(size(labelvolume), 'uint8');

num = length(label);

for i = 1:num
    label3D(find(labelvolume==label(i))) = 1;
end
return;