function write_color_lut(filename, map);
% function write_color_lut(filename, map);

n_lut = size(map, 1);

% interpolate to 4096 slots
p = map;
map = interp1(0:n_lut-1, p, 0:((n_lut-1)/4095):n_lut-1, 'pchip');

% scale make map between 0-255
map = 255*map/max(map(:));
map = uint8(map);
red = map(:,1);
green = map(:,2);
blue = map(:,3);

fid=fopen(filename,'w');
fprintf(fid,'%s\n','Siemens AG, Copyright 2008');
fprintf(fid,'%s\n','Numaris 4, VD13A');
fprintf(fid,'%s\n','BoldColorScale');
fprintf(fid,'%s\n','Version : 1.0');
fprintf(fid,'%s\n',num2str(size(map,1)));
for i = 1:size(map,1)
    fprintf(fid,'%s%s%s%s%s%s%s%s\n','[',num2str(i-1),'] ',num2str(red(i)),' , ',num2str(green(i)),' , ',num2str(blue(i)));    
end
fclose(fid); 

return
