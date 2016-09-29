
function [xi, yi, zi] = generateSampledCoordinatesPerfusionPhantomSAX(data)

[rows,cols,slices]=size(data);

data1 = data;
x = 1:rows;
y = 1:cols;
z = 1:slices;
offset=[1:slices];

for k=1:length(offset)
    k
    tic
    hsp = surf(1:rows, 1:cols, zeros(rows,cols)+offset(k));
    rotate(hsp,[0,1,0],-20)
    rotate(hsp,[0,0,1],-40)
%     rotate(hsp,[0,1,0],-20)
%     rotate(hsp,[0,1,0],-90)
    xi(:,:,k) = get(hsp,'XData');
    yi(:,:,k) = get(hsp,'YData');
    zi(:,:,k) = get(hsp,'ZData');
    close
    data1(:,:,k)=interp3(x,y,z,data,xi(:,:,k),yi(:,:,k),zi(:,:,k), 'cubic');
    toc
end
shg;

%%
% [rows,cols,slices]=size(data1);
% % [x,y,z] = meshgrid(1:rows,1:cols,1:slices);
% x = 1:rows;
% y = 1:cols;
% z = 1:slices;
% clear data2
% for k=1:length(offsetX)
%     k
%     hsp = surf([1:rows],1:cols,zeros(rows,cols));
%     rotate(hsp,[0,1,0],90)
%     xi2(:,:,k) = get(hsp,'XData');
%     yi2(:,:,k) = get(hsp,'YData');
%     zi2(:,:,k) = get(hsp,'ZData');
%     xi2(:,:,k) = xi2(:,:,k) + offsetX(k);
%     close
% end
% shg;

xi = single(xi);
yi = single(yi);
zi = single(zi);
% xi2 = single(xi2);
% yi2 = single(yi2);
% zi2 = single(zi2);
