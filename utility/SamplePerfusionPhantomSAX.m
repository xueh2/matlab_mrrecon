
function sax = SamplePerfusionPhantomSAX(data, offsetX, xi, yi, zi)

[rows,cols,slices]=size(data);

data1 = data;
x = 1:rows;
y = 1:cols;
z = 1:slices;
offset=[1:slices];

data1 = interp3(x,y,z, double(data), xi(:),yi(:),zi(:), 'cubic');
data1 = reshape(data1, size(data));

% for k=1:length(offset)
%     k
%     tic
% %     hsp = surf(1:rows, 1:cols, zeros(rows,cols)+offset(k));
% %     rotate(hsp,[0,1,0],-20)
% %     rotate(hsp,[0,0,1],-40)
% %     xi = get(hsp,'XData');
% %     yi = get(hsp,'YData');
% %     zi = get(hsp,'ZData');
% %     close
%     data1(:,:,k)=interp3(x,y,z,data,xi(:,:,k),yi(:,:,k),zi(:,:,k), 'cubic');
%     toc
% end
% shg;

%%
% [rows,cols,slices]=size(data1);
% % [x,y,z] = meshgrid(1:rows,1:cols,1:slices);
% x = 1:rows;
% y = 1:cols;
% z = 1:slices;
% 
% data2 = interp3(x,y,z,data1,xi2(:),yi2(:),zi2(:), 'cubic');
% data2 = reshape(data2, size(xi2));

data2 = data1(:, offsetX(1):offsetX(2), :);
data3 = permute(data2, [1 3 2]);
data3 = permute(data3, [2 1 3]);
data3 = flipdim(data3, 1);

% for k=1:length(offsetX)
%     k
% %     hsp = surf([1:rows],1:cols,zeros(rows,cols));
% %     rotate(hsp,[0,1,0],90)
% %     xi2 = get(hsp,'XData');
% %     yi2 = get(hsp,'YData');
% %     zi2 = get(hsp,'ZData');
% %     xi2 = xi2 + offsetX(k);
% %     close
%     data2(:,:,k)=interp3(x,y,z,data1,xi2(:,:,k),yi2(:,:,k),zi2(:,:,k), 'cubic');
% end
% shg;

data4 = data3;
mark = isnan(data3);
data4(find(mark>0)) = 0;

% data3 = permute(data2, [2 1 3]);
% data3 = flipdim(data3, 2);
sax = single(data4);

