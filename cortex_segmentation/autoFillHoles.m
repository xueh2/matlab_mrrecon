
function data_noholes = autoFillHoles(data, header)
% fill the holes within the binary volume

data_noholes = data;

pp = data;
pp = 1 - data;

label = 1;
[l,num] = bwlabeln(pp,6);

xsize=header.xsize;
ysize=header.ysize;
zsize=header.zsize;
volumes=zeros(num,1);
total = xsize*ysize*zsize;
for i=1:total
    if (l(i)>0)
        volumes(l(i))=volumes(l(i))+1;
    end
end
vd = sort(volumes, 'descend');
[largest,index] = max(volumes);

data_noholes(find(l~=index(1))) = 1;

return