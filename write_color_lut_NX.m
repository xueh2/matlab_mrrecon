function write_color_lut_NX(filename, map, num_entries);
% function write_color_lut_NX(filename, map);

n_lut = size(map, 1);

if(nargin<3)
    num_entries = 4096;
end

% interpolate to 4096 slots
% p = map;
% map = interp1(0:255, p, 0:(255/(num_entries-1)):255, 'PCHIP');

p = map;
map = interp1(0:n_lut-1, p, 0:((n_lut-1)/(num_entries-1)):n_lut-1, 'pchip');

% scale make map between 0-255
map = 65535*map/max(map(:));
map = uint16(map);
red = map(:,1);
green = map(:,2);
blue = map(:,3);

fid=fopen(filename,'w');
fprintf(fid,'%s\n','Siemens Healthcare GmbH, Copyright 2015');
fprintf(fid,'%s\n','Numaris/X');
fprintf(fid,'%s\n','PerfusionColorScale');
fprintf(fid,'%s\n','Version : 1.0');
fprintf(fid,'%s\n',num2str(size(map,1)));
for i = 1:size(map,1)
    fprintf(fid,'%s%s%s%s%s%s%s%s\n','[',num2str(i-1),'] ',num2str(red(i)),' , ',num2str(green(i)),' , ',num2str(blue(i)));    
end
fclose(fid); 

return
